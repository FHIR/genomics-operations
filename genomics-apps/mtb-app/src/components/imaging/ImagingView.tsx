import { useEffect, useMemo, useRef, useState } from 'react';
import { ChevronDown, ChevronUp, Maximize2, TriangleAlert } from 'lucide-react';
import { ClinicalState } from '@/lib/useClinicalRecord';
import { useBannerHeight } from '@/lib/useBannerHeight';
import { formatDate } from '@/lib/clinicalSummary';
import { resourceKey } from '@/lib/caseSummary';
import { buildImageCatalog, idcViewerUrl, ImageEntry, reportTitle } from '@/lib/imaging';
import { findDiagnosticSlide } from '@/services/idcService';

// Set when a report's image link in the Case Summary is clicked; a new nonce repeats the jump for the same report
export interface ImagingFocusRequest {
  reportKey: string;
  nonce: number;
}

interface ImagingViewProps {
  subjectId: string;
  clinicalState: ClinicalState;
  onRetry: () => void;
  focusRequest: ImagingFocusRequest | null;
}

// The viewer fills the window below the sticky patient banner, but never gets shorter than this
const MIN_VIEWER_HEIGHT = 600;
const BANNER_GAP = 24;

// Date, type, study, report
const INDEX_COLUMNS = 'grid-cols-[6rem_4.5rem_minmax(0,3fr)_minmax(0,2fr)]';

const TYPE_CHIP_CLASSES = {
  slide: 'bg-pink-100 text-pink-800',
  radiology: 'bg-sky-100 text-sky-800',
};

function TypeChip({ entry }: { entry: ImageEntry }) {
  return (
    <span className={`whitespace-nowrap rounded-full px-2 py-px text-xs font-semibold ${entry.isSlide ? TYPE_CHIP_CLASSES.slide : TYPE_CHIP_CLASSES.radiology}`}>
      {entry.typeLabel}
    </span>
  );
}

function studyLabel(entry: ImageEntry) {
  return entry.isSlide && entry.slideCount > 1 ? `${entry.description} · ${entry.slideCount} slides` : entry.description;
}

function StatusMessage({ subjectId, clinicalState, onRetry }: Omit<ImagingViewProps, 'focusRequest'>) {
  return (
    <div className="rounded-2xl border border-slate-200 bg-white p-8 text-center text-slate-600 shadow-sm">
      {!subjectId.trim() && 'Choose a patient in the banner to see their images.'}
      {subjectId.trim() && clinicalState.status === 'loading' && 'Loading clinical data…'}
      {clinicalState.status === 'none' && `No clinical data found for MRN ${subjectId}, so there are no images to show.`}
      {clinicalState.status === 'error' && (
        <span className="inline-flex flex-col items-center gap-3">
          <span className="text-red-700" title={clinicalState.message}>Clinical data could not be loaded.</span>
          <button type="button" onClick={onRetry} className="rounded-md border border-slate-300 px-3 py-1.5 text-sm font-medium hover:bg-slate-50">Retry</button>
        </span>
      )}
    </div>
  );
}

export default function ImagingView({ subjectId, clinicalState, onRetry, focusRequest }: ImagingViewProps) {
  const record = clinicalState.status === 'loaded' ? clinicalState.record : null;
  const catalog = useMemo(() => (record ? buildImageCatalog(record.resources) : []), [record]);

  const [openKey, setOpenKey] = useState<string | null>(null);
  const [isExpanded, setIsExpanded] = useState(true);
  // The report whose images were asked for from the Case Summary; its rows stay highlighted
  const [highlightedReportKey, setHighlightedReportKey] = useState<string | null>(null);
  // The viewer URL once known; slide studies first look up their FFPE diagnostic slide
  const [viewer, setViewer] = useState<{ key: string; url?: string } | null>(null);
  const catalogRef = useRef<HTMLElement>(null);
  const frameBoxRef = useRef<HTMLDivElement>(null);
  const bannerHeight = useBannerHeight();

  const openEntry = catalog.find((entry) => entry.key === openKey);
  // The table can only collapse while an image is open
  const showTable = isExpanded || !openEntry;

  const open = (key: string) => {
    setOpenKey(key);
    setIsExpanded(false);
    // Once the table has collapsed, bring its one-line bar up under the banner with the viewer right below it
    window.requestAnimationFrame(() => catalogRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' }));
  };

  // A report's image link: expand the table, highlight that report's images and open its newest one
  useEffect(() => {
    if (!focusRequest || !record) {
      return;
    }
    const entries = catalog.filter((entry) => entry.reports.some((report) => resourceKey(report) === focusRequest.reportKey));
    setHighlightedReportKey(focusRequest.reportKey);
    setOpenKey(entries[0]?.key ?? null);
    setIsExpanded(true);
    // The link sits at the bottom of the Case Summary; start the Imaging tab at its table
    window.requestAnimationFrame(() => catalogRef.current?.scrollIntoView({ block: 'start' }));
    // Only a new request should move the view, not a reload of the same record
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [focusRequest?.nonce, record]);

  useEffect(() => {
    if (!openEntry) {
      setViewer(null);
      return;
    }
    if (!openEntry.isSlide || !openEntry.dicomWebBase || !openEntry.studyUid) {
      setViewer({ key: openEntry.key, url: idcViewerUrl(openEntry) });
      return;
    }
    let isCurrent = true;
    setViewer(null);
    findDiagnosticSlide(openEntry.dicomWebBase, openEntry.studyUid)
      .catch((error) => {
        // Fall back to the slide IDC's viewer picks
        console.warn(`Could not look up the diagnostic slide for study ${openEntry.studyUid}:`, error);
        return undefined;
      })
      .then((seriesUid) => {
        if (isCurrent) {
          setViewer({ key: openEntry.key, url: idcViewerUrl(openEntry, seriesUid) });
        }
      });
    return () => {
      isCurrent = false;
    };
  }, [openEntry]);

  if (!record) {
    return <StatusMessage subjectId={subjectId} clinicalState={clinicalState} onRetry={onRetry} />;
  }

  if (catalog.length === 0) {
    return (
      <div className="rounded-2xl border border-slate-200 bg-white p-8 text-center text-slate-600 shadow-sm">
        No images in this patient&apos;s record.
      </div>
    );
  }

  const viewerUrl = viewer && viewer.key === openEntry?.key ? viewer.url : undefined;
  const isResolving = Boolean(openEntry) && (!viewer || viewer.key !== openEntry?.key);
  const fromIdc = catalog.some((entry) => entry.dicomWebBase);
  const viewerHeight = `max(${MIN_VIEWER_HEIGHT}px, calc(100vh - ${bannerHeight + BANNER_GAP * 2}px))`;

  return (
    <div className="flex flex-col gap-4">
      <section ref={catalogRef} aria-label="Image catalog" style={{ scrollMarginTop: bannerHeight + BANNER_GAP }} className="rounded-2xl border border-slate-200 bg-white px-4 pb-2.5 pt-2 shadow-[0_1px_3px_rgba(15,23,42,0.06),0_4px_12px_rgba(15,23,42,0.04)]">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <button
            type="button"
            onClick={() => setIsExpanded((current) => !current)}
            disabled={!openEntry}
            aria-expanded={showTable}
            aria-controls="image-index"
            className="inline-flex min-h-10 items-center gap-2.5 py-1 text-left disabled:cursor-default"
          >
            {openEntry && (
              <span className="flex h-7 w-7 shrink-0 items-center justify-center rounded-lg bg-slate-100 ring-1 ring-inset ring-slate-200" aria-hidden="true">
                {showTable ? <ChevronUp size={16} strokeWidth={2.25} /> : <ChevronDown size={16} strokeWidth={2.25} />}
              </span>
            )}
            <span className="text-base font-semibold text-slate-900">Images ({catalog.length})</span>
            {!showTable && openEntry && (
              <span className="text-sm text-slate-600">
                Showing <b className="text-slate-900">{openEntry.typeLabel} · {studyLabel(openEntry)}</b> · {formatDate(openEntry.date) ?? 'no date'}
              </span>
            )}
          </button>
          {fromIdc && (
            <span className="text-xs text-slate-500">
              Images: <a href="https://portal.imaging.datacommons.cancer.gov/" target="_blank" rel="noopener noreferrer" className="text-amber-700 hover:underline">NCI Imaging Data Commons</a> (CC BY)
            </span>
          )}
        </div>

        {showTable && (
          // Each row is one button (as in Lab results): click anywhere on it to open that study; the open row is
          // the darker amber one, and a report's other studies stay in light amber
          <div id="image-index" className="mt-1 overflow-x-auto">
            <div className="min-w-[640px] text-[13px]">
              <div className={`grid ${INDEX_COLUMNS} gap-3 border-b border-slate-200 px-2 py-1 text-xs font-semibold text-slate-500`}>
                <span>Date</span><span>Type</span><span>Study</span><span>Report</span>
              </div>
              {catalog.map((entry) => {
                const isOpen = entry.key === openKey;
                const isHighlighted = Boolean(highlightedReportKey) && entry.reports.some((report) => resourceKey(report) === highlightedReportKey);
                return (
                  <button
                    key={entry.key}
                    type="button"
                    aria-current={isOpen ? 'true' : undefined}
                    onClick={() => {
                      if (!isHighlighted) {
                        setHighlightedReportKey(null);
                      }
                      open(entry.key);
                    }}
                    className={`grid w-full ${INDEX_COLUMNS} items-center gap-3 border-b border-slate-100 px-2 py-1 text-left ${isOpen
                      ? 'bg-amber-100 shadow-[inset_3px_0_0_#d97706]'
                      : isHighlighted ? 'bg-amber-50 hover:bg-amber-100/60' : 'hover:bg-slate-50'}`}
                  >
                    <span className="whitespace-nowrap text-slate-700">{formatDate(entry.date) ?? '—'}</span>
                    <span><TypeChip entry={entry} /></span>
                    <span className="font-semibold text-slate-900">{studyLabel(entry)}</span>
                    <span className="text-slate-600">{entry.reports.length > 0 ? entry.reports.map(reportTitle).join('; ') : '—'}</span>
                  </button>
                );
              })}
            </div>
          </div>
        )}
      </section>

      {openEntry && (
        <section
          aria-label="Image viewer"
          className="flex flex-col gap-2.5 rounded-2xl border border-slate-200 bg-white p-4 shadow-[0_1px_3px_rgba(15,23,42,0.06),0_4px_12px_rgba(15,23,42,0.04)]"
        >
          <div className="flex flex-wrap items-baseline justify-between gap-2">
            <div className="flex flex-wrap items-baseline gap-x-2.5 gap-y-1">
              <h2 className="text-base font-semibold text-slate-900">{studyLabel(openEntry)}</h2>
              <span className="text-sm text-slate-500">
                {openEntry.typeLabel} · {formatDate(openEntry.date) ?? 'no date'}
                {openEntry.reports.length > 0 && ` · Report: ${openEntry.reports.map(reportTitle).join('; ')}`}
              </span>
            </div>
            {viewerUrl && (
              <button
                type="button"
                onClick={() => frameBoxRef.current?.requestFullscreen?.()}
                className="inline-flex items-center gap-1 text-sm font-semibold text-amber-700 hover:text-amber-900"
              >
                <Maximize2 size={14} aria-hidden="true" />
                Full screen
              </button>
            )}
          </div>

          {openEntry.imageDate && openEntry.imageDate !== openEntry.date && (
            <div role="note" className="flex items-start gap-2.5 rounded-lg border border-amber-500 bg-amber-50 px-3.5 py-2.5 text-sm text-amber-900">
              <TriangleAlert size={18} aria-hidden="true" className="mt-px shrink-0 text-amber-700" />
              <span>
                <b>Report date is synthetically adjusted, and differs from actual image date.</b>
                {' '}Report: {formatDate(openEntry.date) ?? 'no date'} · Image: {formatDate(openEntry.imageDate)}
              </span>
            </div>
          )}

          <div ref={frameBoxRef} className="overflow-hidden rounded-lg border border-slate-300 bg-slate-950" style={{ height: viewerHeight }}>
            {viewerUrl ? (
              <iframe
                key={viewerUrl}
                src={viewerUrl}
                title={`Image viewer: ${studyLabel(openEntry)}`}
                allow="fullscreen"
                className="h-full w-full border-0"
              />
            ) : (
              <div className="flex h-full items-center justify-center p-8 text-center text-sm text-slate-300">
                {isResolving ? 'Opening the viewer…' : 'These images aren’t in the NCI Imaging Data Commons, so there is no viewer for them.'}
              </div>
            )}
          </div>
        </section>
      )}
    </div>
  );
}
