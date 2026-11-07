import { EvidenceLevel, ReportedImplication, ReportedVariant, ReportedVariantsData } from '@/lib/reportedVariants';
import { SectionCard } from './SectionCard';

// The Reported variants card:
// one group per genetic report, each variant with its therapeutic implications and their evidence levels

const EVIDENCE_TONES: Partial<Record<EvidenceLevel, string>> = {
  A: 'bg-emerald-100 text-emerald-800 ring-emerald-200',
  B: 'bg-sky-100 text-sky-800 ring-sky-200',
};
const DEFAULT_TONE = 'bg-slate-100 text-slate-700 ring-slate-200';

function EvidencePill({ implication }: { implication: ReportedImplication }) {
  if (!implication.evidence) {
    return null;
  }
  const hover = [implication.evidence, implication.phenotype].filter(Boolean).join(' · ');
  return (
    <span
      title={hover}
      className={`inline-flex h-6 shrink-0 items-center gap-1 whitespace-nowrap rounded-full px-2 text-xs ring-1 ring-inset ${(implication.level && EVIDENCE_TONES[implication.level]) || DEFAULT_TONE}`}
    >
      <span className="font-medium">Evidence</span>
      {/* A level that doesn't start with a letter is shown as recorded */}
      <span className="font-bold">{implication.level ?? implication.evidence}</span>
    </span>
  );
}

function ImplicationRow({ implication }: { implication: ReportedImplication }) {
  return (
    <div className="flex items-center gap-2.5 rounded-lg bg-slate-50 px-2.5 py-2 ring-1 ring-inset ring-slate-200">
      <div className="flex min-w-0 flex-1 flex-wrap items-baseline gap-x-2 gap-y-0.5">
        {implication.drug && <span className="font-semibold text-slate-900">{implication.drug}</span>}
        {implication.implication && <span className="text-slate-700">{implication.implication}</span>}
      </div>
      <EvidencePill implication={implication} />
    </div>
  );
}

function VariantItem({ variant }: { variant: ReportedVariant }) {
  return (
    <li className="flex flex-col gap-1.5">
      <div className="flex flex-col gap-0.5">
        <div title={variant.details.join(' · ') || undefined} className="flex flex-wrap items-baseline gap-x-2 gap-y-1">
          {variant.gene && <span className="text-[15px] font-bold text-slate-900">{variant.gene}</span>}
          {variant.change && <span className="break-all text-slate-800">{variant.change}</span>}
        </div>
        {variant.origin && <span className="text-xs text-slate-500">{variant.origin}</span>}
      </div>
      {variant.implications.map((implication) => <ImplicationRow key={implication.key} implication={implication} />)}
      {variant.implications.length === 0 && <span className="text-[13px] text-slate-500">No therapeutic implication reported</span>}
    </li>
  );
}

function GroupHeader({ title, date, details }: { title: string; date?: string; details: string[] }) {
  return (
    <div
      title={details.join(' · ') || undefined}
      className="flex flex-wrap items-baseline justify-between gap-x-2 gap-y-0.5 rounded-lg bg-slate-100 px-2.5 py-1.5"
    >
      <span className="text-xs font-semibold uppercase tracking-wide text-slate-600">{title}</span>
      {date && <span className="text-xs font-semibold text-slate-900">{date}</span>}
    </div>
  );
}

interface ReportedVariantsSectionProps {
  data: ReportedVariantsData;
  // Opens the Genomics Explorer tab, which has the patient's full genomic data set
  onOpenGenomicsExplorer: () => void;
  className?: string;
}

export default function ReportedVariantsSection({ data, onOpenGenomicsExplorer, className = '' }: ReportedVariantsSectionProps) {
  const { groups, otherImplications } = data;
  const footer = (
    <button
      type="button"
      onClick={onOpenGenomicsExplorer}
      className="rounded text-xs font-semibold text-emerald-700 hover:text-emerald-900 hover:underline focus:outline-none focus:ring-2 focus:ring-emerald-300"
    >
      See full genomic data in Genomics Explorer →
    </button>
  );
  return (
    <SectionCard section="variants" fillHeight footer={footer} className={className}>
      {groups.length === 0 && otherImplications.length === 0 ? (
        <p className="text-sm text-slate-500">No genetic report in this record.</p>
      ) : (
        <div className="flex flex-col gap-5 text-sm">
          {groups.map((group) => (
            <div key={group.key} className="flex flex-col gap-3">
              <GroupHeader title={group.title} date={group.date} details={group.details} />
              {group.variants.length > 0 ? (
                <ul className="flex flex-col gap-3.5">
                  {group.variants.map((variant) => <VariantItem key={variant.key} variant={variant} />)}
                </ul>
              ) : (
                <div className="flex flex-col gap-1.5">
                  {group.conclusion && <p className="leading-relaxed text-slate-800">{group.conclusion}</p>}
                  <span className="text-[13px] text-slate-500">This report has no structured variants. See Reports for the full text.</span>
                </div>
              )}
            </div>
          ))}
          {otherImplications.length > 0 && (
            <div className="flex flex-col gap-3">
              <GroupHeader title="Other implications" details={[]} />
              {otherImplications.map((implication) => <ImplicationRow key={implication.key} implication={implication} />)}
            </div>
          )}
        </div>
      )}
    </SectionCard>
  );
}
