import { useEffect, useState } from 'react';
import { Info } from 'lucide-react';

// Search status as a row of chips, one per search term (mockup: Search Workspace redesign, states B–F).
// Terms are searched one at a time in order, so the first term still marked 'searching' is the one in
// progress and the rest are waiting.

// Past this many terms the chips stop after about two rows, with "Show all" to expand
const CHIP_CAP = 12;

interface SearchStatusProps {
  // Search term → 'searching' | 'found' | 'no results' | 'not_on_server' | 'error' | 'cancelled'
  searchStatus: Record<string, string>;
  // The genomics server has no data for this patient, so results come only from the genetic report
  noServerData: boolean;
  invalidRanges: string[];
  // Variants found per search term
  variantCounts: Record<string, number>;
  onCancelSearch: () => void;
}

type ChipKind = 'found' | 'searching' | 'waiting' | 'none' | 'error' | 'cancelled';

interface Chip {
  term: string;
  kind: ChipKind;
  label: string;
}

const CHIP_CLASSES: Record<ChipKind, string> = {
  found: 'bg-emerald-50 text-emerald-800 ring-emerald-200',
  searching: 'bg-blue-50 text-blue-800 ring-blue-200',
  waiting: 'bg-white text-slate-400 ring-slate-200',
  none: 'bg-slate-50 text-slate-500 ring-slate-200',
  error: 'bg-red-50 text-red-800 ring-red-200',
  cancelled: 'bg-orange-50 text-orange-800 ring-orange-200',
};

function buildChips(searchStatus: Record<string, string>, variantCounts: Record<string, number>, noServerData: boolean): Chip[] {
  const fromReport = noServerData ? ' (report)' : '';
  const current = Object.keys(searchStatus).find((term) => searchStatus[term] === 'searching');
  return Object.entries(searchStatus).map(([term, status]) => {
    switch (status) {
      case 'found': {
        const count = variantCounts[term] ?? 0;
        return { term, kind: 'found', label: count > 0 ? `✓ ${term} · ${count} variant${count === 1 ? '' : 's'}${fromReport}` : `✓ ${term} · Results found` };
      }
      case 'searching':
        return term === current
          ? { term, kind: 'searching', label: `◌ ${term} · Searching…` }
          : { term, kind: 'waiting', label: `${term} · Waiting` };
      case 'no results':
        return { term, kind: 'none', label: `○ ${term} · No results found` };
      case 'not_on_server':
        return { term, kind: 'none', label: `${term} · Not on genomics server` };
      case 'cancelled':
        return { term, kind: 'cancelled', label: `${term} · Cancelled` };
      default:
        return { term, kind: 'error', label: `✕ ${term} · Error processing` };
    }
  });
}

/**
 * One-line summary of a search, e.g. "46 of 83 searched · 3 found · 43 no results found", or with no data on the
 * genomics server "11 searched · 1 found in the report · 10 not on the genomics server".
 */
export function searchSummaryText(searchStatus: Record<string, string>, noServerData = false) {
  const statuses = Object.values(searchStatus);
  const count = (...values: string[]) => statuses.filter((status) => values.includes(status)).length;
  const done = statuses.length - count('searching');
  const parts = [
    count('searching') > 0 ? `${done} of ${statuses.length} searched` : `${statuses.length} searched`,
    noServerData ? `${count('found')} found in the report` : `${count('found')} found`,
  ];
  if (!noServerData || count('no results') > 0) {
    parts.push(`${count('no results')} no results found`);
  }
  if (count('not_on_server') > 0) {
    parts.push(`${count('not_on_server')} not on the genomics server`);
  }
  if (count('error') > 0) {
    parts.push(`${count('error')} with errors`);
  }
  if (count('cancelled') > 0) {
    parts.push(`${count('cancelled')} cancelled`);
  }
  return parts.join(' · ');
}

function ChipView({ chip }: { chip: Chip }) {
  return (
    <span className={`whitespace-nowrap rounded-full px-2.5 py-1 text-xs ring-1 ring-inset ${CHIP_CLASSES[chip.kind]}`}>{chip.label}</span>
  );
}

export default function SearchStatus({ searchStatus, noServerData, invalidRanges, variantCounts, onCancelSearch }: SearchStatusProps) {
  const [showAll, setShowAll] = useState(false);
  const terms = Object.keys(searchStatus);
  const termKey = terms.join('\n');

  // A new search starts collapsed again
  useEffect(() => {
    setShowAll(false);
  }, [termKey]);

  if (terms.length === 0 && invalidRanges.length === 0) {
    return null;
  }

  const chips = buildChips(searchStatus, variantCounts, noServerData);
  const isSearching = chips.some((chip) => chip.kind === 'searching' || chip.kind === 'waiting');
  const doneCount = chips.filter((chip) => chip.kind !== 'searching' && chip.kind !== 'waiting').length;
  const isLong = chips.length > CHIP_CAP;

  // Capped view: results, the term in progress and errors always show; then the most recent "No results found"
  const priority = chips.filter((chip) => ['found', 'searching', 'error'].includes(chip.kind));
  const recent = chips.filter((chip) => chip.kind === 'none' || chip.kind === 'cancelled').reverse();
  const capped = [...priority, ...recent.slice(0, Math.max(0, CHIP_CAP - priority.length))];
  const hiddenCount = chips.length - capped.length;

  return (
    <div className="flex flex-col gap-2.5 border-t border-emerald-100 pt-3">
      {noServerData && (
        <div className="flex items-start gap-2 rounded-lg bg-amber-50 px-3 py-2.5 text-sm text-amber-900 ring-1 ring-inset ring-amber-300">
          <Info size={16} className="mt-0.5 flex-shrink-0" aria-hidden="true" />
          <span>
            <span className="font-semibold">The genomics server has no data for this patient.</span>{' '}
            Results below come only from the patient&apos;s genetic report.
          </span>
        </div>
      )}
      <div className="flex flex-wrap items-baseline gap-x-2 gap-y-1 text-sm">
        <span className="font-semibold text-slate-700">Status</span>
        <span className="text-slate-600">{searchSummaryText(searchStatus, noServerData)}</span>
        <span className="flex-1" />
        {isSearching && (
          <button type="button" onClick={onCancelSearch} className="font-semibold text-red-700 hover:underline">
            Cancel search
          </button>
        )}
      </div>

      {isSearching && chips.length > 1 && (
        <div className="h-1.5 overflow-hidden rounded-full bg-slate-100" role="progressbar" aria-valuemin={0} aria-valuemax={chips.length} aria-valuenow={doneCount}>
          <div className="h-full rounded-full bg-emerald-500 transition-[width] duration-500" style={{ width: `${(doneCount / chips.length) * 100}%` }} />
        </div>
      )}

      {isLong && showAll ? (
        <>
          {/* Every term in search order, in a box of fixed height that scrolls */}
          <div className="flex max-h-36 flex-wrap content-start gap-2 overflow-y-auto rounded-lg border border-emerald-100 bg-white p-2.5">
            {chips.map((chip) => <ChipView key={chip.term} chip={chip} />)}
          </div>
          <button type="button" onClick={() => setShowAll(false)} className="self-start text-sm font-semibold text-blue-600 hover:underline">
            Show fewer
          </button>
        </>
      ) : (
        <div className="flex flex-wrap items-center gap-2">
          {(isLong ? capped : chips).map((chip) => <ChipView key={chip.term} chip={chip} />)}
          {isLong && hiddenCount > 0 && (
            <button type="button" onClick={() => setShowAll(true)} className="ml-1 text-sm font-semibold text-blue-600 hover:underline">
              +{hiddenCount} more · Show all {chips.length}
            </button>
          )}
        </div>
      )}

      {invalidRanges.length > 0 && (
        <p className="text-sm text-red-700">
          Unprocessable: {invalidRanges.join(', ')}. Unrecognized gene or malformed chromosomal range.
        </p>
      )}
    </div>
  );
}
