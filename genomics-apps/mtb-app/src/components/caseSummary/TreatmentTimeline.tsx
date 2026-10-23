import { CSSProperties, Fragment, useEffect, useRef, useState } from 'react';
import { MarkerTone, TimelineBar, TimelineMarker, TreatmentTimeline as Timeline } from '@/lib/treatmentTimeline';
import { SectionCard } from './SectionCard';

// Lane label column width; the plot takes the rest
const LABEL_COLUMN = '15rem';
// Past this share of the axis, labels go to the left of their mark so they don't run off the edge
const RIGHT_EDGE = 70;
// Bars at least this wide carry their label inside
const LABEL_INSIDE_WIDTH = 30;
// Marker rows: labels that would overlap drop to another line within the lane
const ROW_HEIGHT = 22;
const LANE_HEIGHT = 32;
const LABEL_GAP = 8;
// How far a marker's box starts before (or ends after) its date: half a chip or half a ◆
const CHIP_OFFSET = 8;
const SYMBOL_OFFSET = 5.6;

const BAR_CLASSES: Record<TimelineBar['kind'], string> = {
  medication: 'bg-blue-600',
  radiotherapy: 'bg-orange-500',
};

const MARKER_CLASSES: Record<MarkerTone, string> = {
  improving: 'bg-emerald-100 text-emerald-900 border-emerald-300',
  stable: 'bg-slate-100 text-slate-800 border-slate-300',
  progression: 'bg-orange-100 text-orange-900 border-orange-300',
  other: 'bg-amber-50 text-amber-900 border-amber-300',
  procedure: 'text-slate-700',
  diagnosis: 'text-violet-800',
};

function position(date: string, timeline: Timeline) {
  const start = Date.parse(timeline.start);
  const span = Date.parse(timeline.end) - start || 1;
  return Math.min(100, Math.max(0, ((Date.parse(date.slice(0, 10)) - start) / span) * 100));
}

// Year ticks between the axis ends, plus the start month
function axisTicks(timeline: Timeline) {
  const ticks: { left: number; label: string }[] = [{
    left: 0,
    label: new Date(`${timeline.start}T00:00:00`).toLocaleDateString('en-US', { month: 'short', year: 'numeric' }),
  }];
  const firstYear = Number(timeline.start.slice(0, 4)) + 1;
  const lastYear = Number(timeline.end.slice(0, 4));
  for (let year = firstYear; year <= lastYear; year++) {
    const left = position(`${year}-01-01`, timeline);
    // Skip a year tick that would collide with the start label
    if (left > 8) {
      ticks.push({ left, label: String(year) });
    }
  }
  return ticks;
}

function Bar({ bar, timeline }: { bar: TimelineBar; timeline: Timeline }) {
  const left = position(bar.start, timeline);
  const width = Math.max(position(bar.end, timeline) - left, 0.6);
  const end = left + width;
  const inside = width >= LABEL_INSIDE_WIDTH;
  const labelStyle: CSSProperties = inside
    ? { left: `${left}%`, width: `${width}%` }
    : end > RIGHT_EDGE
      ? { right: `${100 - left}%` }
      : { left: `${end}%` };
  return (
    <>
      <div
        title={bar.tooltip}
        className={`absolute top-1/2 h-3.5 -translate-y-1/2 rounded-sm ${BAR_CLASSES[bar.kind]}`}
        style={{ left: `${left}%`, width: `${width}%` }}
      />
      <span
        title={bar.tooltip}
        className={`absolute top-1/2 -translate-y-1/2 whitespace-nowrap px-1.5 text-xs ${inside ? 'truncate font-medium text-white' : 'text-slate-600'}`}
        style={labelStyle}
      >
        {bar.label}
      </span>
    </>
  );
}

let measureContext: CanvasRenderingContext2D | null = null;

function textWidth(text: string, font: string) {
  measureContext ??= document.createElement('canvas').getContext('2d');
  if (!measureContext) {
    return text.length * 6.5;
  }
  measureContext.font = font;
  return measureContext.measureText(text).width;
}

/**
 * The line (0, 1, …) each marker sits on within its lane. Labels are measured at the plot's width, and
 * one that would overlap an earlier label goes to the first line where it fits. A marker without a label
 * shares the line of the labeled marker before it (its group).
 */
function markerRows(markers: TimelineMarker[], timeline: Timeline, plotWidth: number, fontFamily: string): number[] {
  if (plotWidth === 0) {
    return markers.map(() => 0);
  }
  const occupied: [number, number][][] = [];
  let groupRow = 0;
  return markers.map((marker) => {
    if (!marker.label) {
      return groupRow;
    }
    const isChip = !marker.symbol;
    const font = `${marker.emphasized ? 'bold ' : ''}12px ${fontFamily}`;
    // Chip padding and border, or the ◆ and its gap
    const width = textWidth(marker.label, font) + (isChip ? 16 : 18) + (marker.emphasized ? 6 : 0);
    const x = (position(marker.date, timeline) / 100) * plotWidth;
    const offset = isChip ? CHIP_OFFSET : SYMBOL_OFFSET;
    const interval: [number, number] = position(marker.date, timeline) > RIGHT_EDGE
      ? [x + offset - width, x + offset]
      : [x - offset, x - offset + width];
    let row = occupied.findIndex((intervals) =>
      intervals.every(([start, end]) => interval[1] + LABEL_GAP <= start || interval[0] >= end + LABEL_GAP));
    if (row === -1) {
      row = occupied.length;
      occupied.push([]);
    }
    occupied[row].push(interval);
    groupRow = row;
    return row;
  });
}

function Marker({ marker, timeline, row }: { marker: TimelineMarker; timeline: Timeline; row: number }) {
  const left = position(marker.date, timeline);
  const isChip = !marker.symbol;
  // Chips (disease status) center on their date; ◆ markers sit on the date with the label beside them
  const anchor = left > RIGHT_EDGE ? 'right' : 'left';
  const top = LANE_HEIGHT / 2 + row * ROW_HEIGHT;
  return (
    <span
      title={marker.tooltip}
      className={`absolute flex -translate-y-1/2 items-center gap-1 whitespace-nowrap text-xs ${anchor === 'right' ? 'flex-row-reverse' : ''}`}
      style={anchor === 'right'
        ? { top, right: `${100 - left}%`, marginRight: isChip ? -CHIP_OFFSET : -SYMBOL_OFFSET }
        : { top, left: `${left}%`, marginLeft: isChip ? -CHIP_OFFSET : -SYMBOL_OFFSET }}
    >
      {isChip ? (
        <span className={`rounded-full border px-1.5 py-px ${MARKER_CLASSES[marker.tone]} ${marker.emphasized ? 'font-bold ring-2 ring-offset-1 ring-slate-400' : ''}`}>
          {marker.label}
        </span>
      ) : (
        <>
          <span className={`text-sm leading-none ${MARKER_CLASSES[marker.tone]}`}>{marker.symbol}</span>
          {marker.label && <span className="text-slate-700">{marker.label}</span>}
        </>
      )}
    </span>
  );
}

// The dashed Today line, drawn in every row so it runs the full height while scrolling
function TodayLine() {
  return <div className="absolute inset-y-0 right-0 border-r-2 border-dashed border-slate-400" aria-hidden="true" />;
}

export default function TreatmentTimeline({ timeline }: { timeline: Timeline | null }) {
  // The plot's width in pixels, for laying out marker labels; follows window resizes
  const plotRef = useRef<HTMLDivElement | null>(null);
  const [plotWidth, setPlotWidth] = useState(0);
  const [fontFamily, setFontFamily] = useState('sans-serif');
  useEffect(() => {
    const element = plotRef.current;
    if (!element) {
      return;
    }
    setFontFamily(getComputedStyle(element).fontFamily);
    const observer = new ResizeObserver(([entry]) => setPlotWidth(entry.contentRect.width));
    observer.observe(element);
    return () => observer.disconnect();
  }, [timeline]);

  if (!timeline) {
    return (
      <SectionCard section="timeline">
        <p className="text-sm text-slate-500">No dated medications, procedures, disease status or diagnoses recorded</p>
      </SectionCard>
    );
  }

  const ticks = axisTicks(timeline);
  const legend = (
    <div className="flex items-center gap-4 text-xs text-slate-600">
      <span className="inline-flex items-center gap-1.5"><span className="h-2.5 w-5 rounded-sm bg-blue-600" />Medication</span>
      <span className="inline-flex items-center gap-1.5"><span className="h-2.5 w-5 rounded-sm bg-orange-500" />Radiotherapy</span>
      <span className="inline-flex items-center gap-1.5"><span className="h-3 border-r-2 border-dashed border-slate-400" />Today</span>
    </div>
  );

  return (
    <SectionCard section="timeline" actions={legend}>
      {/* Fixed maximum height (about 10 lanes); the axis header stays visible while scrolling */}
      {/* Scrolls sideways on narrow windows, so the time axis keeps a readable width */}
      <div className="max-h-[24rem] overflow-auto rounded-lg border border-slate-200">
        <div className="min-w-[48rem]">
          <div className="sticky top-0 z-10 grid border-b border-slate-200 bg-slate-50" style={{ gridTemplateColumns: `${LABEL_COLUMN} 1fr` }}>
            <div />
            <div ref={plotRef} className="relative mr-24 h-7 text-xs text-slate-500">
              {ticks.map((tick) => (
                <span key={tick.label} className="absolute top-1.5 border-l border-slate-300 pl-1" style={{ left: `${tick.left}%` }}>{tick.label}</span>
              ))}
              <span className="absolute top-1.5 right-0 translate-x-full pl-1 font-semibold text-slate-700">Today</span>
              <TodayLine />
            </div>
          </div>
          {timeline.lanes.map((lane, index) => {
            const showHeading = index === 0 || timeline.lanes[index - 1].group !== lane.group;
            const rows = markerRows(lane.markers, timeline, plotWidth, fontFamily);
            const height = LANE_HEIGHT + Math.max(0, ...rows) * ROW_HEIGHT;
            return (
              <Fragment key={lane.key}>
                {showHeading && (
                  // A shaded band marks the start of each block of lanes
                  <div className={`grid border-b border-slate-200 bg-slate-100 ${index > 0 ? 'border-t' : ''}`} style={{ gridTemplateColumns: `${LABEL_COLUMN} 1fr` }}>
                    <div className="px-3 py-1 text-xs font-bold uppercase tracking-wide text-slate-700">{lane.group}</div>
                    <div className="relative mr-24"><TodayLine /></div>
                  </div>
                )}
                <div className="grid hover:bg-slate-50" style={{ gridTemplateColumns: `${LABEL_COLUMN} 1fr` }}>
                  <div className="px-3 py-1 text-sm leading-tight text-slate-900">{lane.label}</div>
                  <div className="relative mr-24" style={{ minHeight: height }}>
                    {ticks.map((tick) => (
                      <div key={tick.label} className="absolute inset-y-0 border-l border-slate-100" style={{ left: `${tick.left}%` }} aria-hidden="true" />
                    ))}
                    <TodayLine />
                    {lane.bars.map((bar) => <Bar key={`${bar.start}-${bar.end}-${bar.label}`} bar={bar} timeline={timeline} />)}
                    {lane.markers.map((marker, markerIndex) => (
                      <Marker key={`${marker.date}-${markerIndex}`} marker={marker} timeline={timeline} row={rows[markerIndex]} />
                    ))}
                  </div>
                </div>
              </Fragment>
            );
          })}
        </div>
      </div>
    </SectionCard>
  );
}
