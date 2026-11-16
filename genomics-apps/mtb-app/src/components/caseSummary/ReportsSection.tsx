import { useState } from 'react';
import { Images } from 'lucide-react';
import { DiagnosticReport } from '@/types/fhir';
import { conceptLabel, formatDate } from '@/lib/clinicalSummary';
import { ReportGroup, reportDate, resourceKey } from '@/lib/caseSummary';
import { SectionCard } from './SectionCard';

// A report's image link: its label ("4 image studies", "3 slides") keyed by the report's resource key, and what a
// click does
export interface ReportImageLinks {
  labels: Map<string, string>;
  onOpen: (reportKey: string) => void;
}

const REPORTS_SHOWN = 3;

// The report's own words: conclusion, else coded conclusions, else a plain-text attachment. Never parsed.
function reportText(report: DiagnosticReport) {
  if (report.conclusion) {
    return report.conclusion;
  }
  const codes = report.conclusionCode?.map((code) => conceptLabel(code)).filter(Boolean);
  if (codes?.length) {
    return codes.join('; ');
  }
  const attachment = report.presentedForm?.find((form) => form.contentType?.startsWith('text/plain') && form.data);
  if (attachment?.data) {
    try {
      return new TextDecoder().decode(Uint8Array.from(atob(attachment.data), (char) => char.charCodeAt(0)));
    } catch {
      return undefined;
    }
  }
  return undefined;
}

function ReportItem({ report, imageLinks }: { report: DiagnosticReport; imageLinks: ReportImageLinks }) {
  const [isExpanded, setIsExpanded] = useState(false);
  const text = reportText(report);
  const key = resourceKey(report);
  const imagesLabel = imageLinks.labels.get(key);
  return (
    <div className="border-t border-slate-200 py-2 first:border-t-0">
      <div className="flex flex-wrap items-center gap-x-2 gap-y-1 text-sm">
        <span>
          <b className="text-slate-900">{conceptLabel(report.code) ?? 'Report'}</b>
          <span className="text-slate-500"> · {formatDate(reportDate(report)) ?? 'no date'}</span>
          {report.status && !['final', 'amended', 'corrected'].includes(report.status) && <span className="text-amber-700"> · {report.status}</span>}
        </span>
        {/* Amber, the Imaging tab's color */}
        {imagesLabel && (
          <button
            type="button"
            onClick={() => imageLinks.onOpen(key)}
            title="Open these images in the Imaging tab"
            className="inline-flex items-center gap-1.5 rounded-full bg-amber-50 px-2.5 py-0.5 text-xs font-semibold text-amber-800 ring-1 ring-inset ring-amber-300 hover:bg-amber-100"
          >
            <Images size={13} aria-hidden="true" />
            {imagesLabel}
          </button>
        )}
      </div>
      {text && (
        <button
          type="button"
          onClick={() => setIsExpanded((current) => !current)}
          aria-expanded={isExpanded}
          title={isExpanded ? 'Show less' : 'Show the full text'}
          className={`mt-0.5 block w-full text-left text-sm text-slate-700 ${isExpanded ? 'whitespace-pre-line' : 'line-clamp-2'}`}
        >
          {text}
        </button>
      )}
    </div>
  );
}

function ReportGroupColumn({ group, imageLinks }: { group: ReportGroup; imageLinks: ReportImageLinks }) {
  const [showAll, setShowAll] = useState(false);
  const shown = showAll ? group.reports : group.reports.slice(0, REPORTS_SHOWN);
  const hiddenCount = group.reports.length - shown.length;
  return (
    <div className="min-w-0 rounded-lg border border-slate-200 p-3">
      <h3 className="mb-1 text-sm font-semibold text-slate-600">{group.label} ({group.reports.length})</h3>
      {shown.map((report) => <ReportItem key={report.id} report={report} imageLinks={imageLinks} />)}
      {(hiddenCount > 0 || showAll) && group.reports.length > REPORTS_SHOWN && (
        <button type="button" onClick={() => setShowAll((current) => !current)} className="mt-1 text-xs font-medium text-blue-600 hover:underline">
          {showAll ? 'Show fewer' : `Show ${hiddenCount} more`}
        </button>
      )}
    </div>
  );
}

export default function ReportsSection({ groups, imageLinks }: { groups: ReportGroup[]; imageLinks: ReportImageLinks }) {
  return (
    <SectionCard section="reports">
      {groups.length === 0 ? <p className="text-sm text-slate-500">None recorded</p> : (
        <div className="grid grid-cols-1 gap-4 lg:grid-cols-2">
          {groups.map((group) => <ReportGroupColumn key={group.label} group={group} imageLinks={imageLinks} />)}
        </div>
      )}
    </SectionCard>
  );
}
