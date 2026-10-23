import { useMemo, useState } from 'react';
import { FhirResource, Observation } from '@/types/fhir';
import { effectiveDate, formatDate } from '@/lib/clinicalSummary';
import {
  interpretationFlag,
  numericValue,
  ObservationGroup,
  observationValueText,
  referenceRangeText,
} from '@/lib/caseSummary';
import { SectionCard } from './SectionCard';

// Trend lines need this many numeric values in one unit
const MIN_TREND_POINTS = 3;
// Shared by the header and the rows; the Test column is widest so long test names fit on one line
const COLUMNS = 'grid-cols-[minmax(0,3.2fr)_minmax(0,1.4fr)_6rem_minmax(0,0.9fr)_minmax(0,0.9fr)_7.5rem]';

interface TrendPoint {
  date: string;
  value: number;
}

function trendPoints(observations: Observation[]): TrendPoint[] | null {
  const points = observations
    .map((observation) => ({ date: effectiveDate(observation), numeric: numericValue(observation) }))
    .filter((point): point is { date: string; numeric: { value: number; unit: string } } => Boolean(point.date && point.numeric));
  const units = new Set(points.map((point) => point.numeric.unit));
  // Results in different units aren't charted (design 4.3, principle 7)
  if (points.length < MIN_TREND_POINTS || units.size > 1) {
    return null;
  }
  return points.map((point) => ({ date: point.date, value: point.numeric.value })).reverse();
}

function scalePoints(points: TrendPoint[], width: number, height: number, padding: number) {
  const days = points.map((point) => Date.parse(point.date));
  const first = Math.min(...days);
  const last = Math.max(...days);
  const values = points.map((point) => point.value);
  const min = Math.min(0, ...values);
  const max = Math.max(...values) || 1;
  return points.map((point, index) => ({
    ...point,
    x: padding + (last > first ? (days[index] - first) / (last - first) : 0.5) * (width - 2 * padding),
    y: height - padding - ((point.value - min) / (max - min)) * (height - 2 * padding),
  }));
}

function Sparkline({ points }: { points: TrendPoint[] }) {
  const scaled = scalePoints(points, 56, 20, 2);
  return (
    <svg width="56" height="20" viewBox="0 0 56 20" aria-hidden="true" className="shrink-0">
      <polyline points={scaled.map((point) => `${point.x},${point.y}`).join(' ')} fill="none" stroke="#2563eb" strokeWidth="1.5" />
    </svg>
  );
}

function TrendChart({ label, points }: { label: string; points: TrendPoint[] }) {
  const width = 640;
  const height = 150;
  const scaled = scalePoints(points, width, height, 26);
  return (
    <svg
      width={width}
      height={height}
      viewBox={`0 0 ${width} ${height}`}
      role="img"
      aria-label={`${label}: ${points.map((point) => `${formatDate(point.date)} ${point.value}`).join(', ')}`}
      className="max-w-full"
    >
      <line x1="26" y1={height - 26} x2={width - 26} y2={height - 26} stroke="#cbd5e1" />
      <polyline points={scaled.map((point) => `${point.x},${point.y}`).join(' ')} fill="none" stroke="#2563eb" strokeWidth="2" />
      {scaled.map((point) => (
        <g key={`${point.date}-${point.value}`}>
          <circle cx={point.x} cy={point.y} r="3.5" fill="#2563eb" />
          <text x={point.x} y={point.y - 9} textAnchor="middle" fontSize="11" fill="#0f172a">{point.value}</text>
          <text x={point.x} y={height - 10} textAnchor="middle" fontSize="10" fill="#64748b">{formatDate(point.date)}</text>
        </g>
      ))}
    </svg>
  );
}

interface LabsSectionProps {
  labs: ObservationGroup[];
  index: Map<string, FhirResource>;
}

export default function LabsSection({ labs, index }: LabsSectionProps) {
  const [sortBy, setSortBy] = useState<'name' | 'date'>('name');
  const [openKey, setOpenKey] = useState<string | null>(null);

  const sorted = useMemo(() => {
    const latestDate = (group: ObservationGroup) => effectiveDate(group.observations[0]) ?? '';
    return [...labs].sort((left, right) => (sortBy === 'date'
      ? latestDate(right).localeCompare(latestDate(left)) || left.label.localeCompare(right.label)
      : left.label.localeCompare(right.label)));
  }, [labs, sortBy]);

  const sortButton = (value: 'name' | 'date', label: string) => (
    <button
      type="button"
      aria-pressed={sortBy === value}
      onClick={() => setSortBy(value)}
      className={`rounded-md border px-3 py-1 text-xs font-semibold ${sortBy === value ? 'border-blue-700 bg-blue-600 text-white' : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'}`}
    >
      {label}
    </button>
  );

  return (
    <SectionCard
      title="Lab & pathology results"
      actions={labs.length > 1 && (
        <div role="group" aria-label="Sort lab and pathology results" className="flex items-center gap-2 text-xs text-slate-500">
          Sort by {sortButton('name', 'Test name')} {sortButton('date', 'Date')}
        </div>
      )}
    >
      {labs.length === 0 ? <p className="text-sm text-slate-500">None recorded</p> : (
        // Fixed maximum height (about 12 rows) with a scroll bar; the header stays visible while scrolling
        <div className="max-h-[26rem] overflow-y-auto rounded-lg border border-slate-200">
          <div className={`sticky top-0 z-10 grid ${COLUMNS} gap-3 bg-slate-50 px-3 py-1.5 text-xs font-semibold uppercase tracking-wide text-slate-500`}>
            <span>Test</span><span>Latest</span><span>Date</span><span>Reference</span><span>Flag</span><span>History</span>
          </div>
          {sorted.map((group) => {
            const latest = group.observations[0];
            const flag = interpretationFlag(latest);
            const points = trendPoints(group.observations);
            const isOpen = openKey === group.key;
            const isPreliminary = latest.status && !['final', 'amended', 'corrected'].includes(latest.status);
            return (
              <div key={group.key} className="border-t border-slate-200">
                <button
                  type="button"
                  aria-expanded={isOpen}
                  onClick={() => setOpenKey(isOpen ? null : group.key)}
                  className={`grid w-full ${COLUMNS} items-center gap-3 px-3 py-1 text-left text-sm leading-snug hover:bg-blue-50/60`}
                >
                  <span title={group.label} className="truncate font-medium text-slate-900">{group.label}</span>
                  <span className="text-slate-800">
                    {observationValueText(latest, index) ?? '—'}
                    {isPreliminary && <span className="ml-1 text-xs text-amber-700">({latest.status})</span>}
                  </span>
                  <span className="text-slate-600">{formatDate(effectiveDate(latest)) ?? '—'}</span>
                  <span className="text-slate-600">{referenceRangeText(latest) ?? ''}</span>
                  <span className={flag?.direction === 'high' ? 'font-semibold text-orange-800' : flag?.direction === 'low' ? 'font-semibold text-blue-800' : 'text-slate-700'}>
                    {flag?.text ?? ''}
                  </span>
                  <span className="flex items-center gap-2 text-xs text-slate-500">
                    {points && <Sparkline points={points} />}
                    {group.observations.length === 1 ? '1 value' : `${group.observations.length} values`}
                  </span>
                </button>
                {isOpen && (
                  <div className="flex flex-wrap items-start gap-6 bg-slate-50 px-4 py-3">
                    {points && <TrendChart label={group.label} points={points} />}
                    <div className="flex flex-col gap-1 text-sm">
                      {group.observations.map((observation) => (
                        <span key={observation.id}>
                          <span className="text-slate-500">{formatDate(effectiveDate(observation)) ?? 'no date'}</span>
                          {' · '}
                          <b className="text-slate-900">{observationValueText(observation, index) ?? '—'}</b>
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            );
          })}
        </div>
      )}
    </SectionCard>
  );
}
