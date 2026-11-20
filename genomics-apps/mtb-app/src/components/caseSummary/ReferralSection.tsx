import { TumorBoardReferral } from '@/lib/caseSummary';
import { CASE_SECTIONS, sectionAnchor } from './sections';
import { SectionCard } from './SectionCard';

// Why the patient is before the tumor board, from the referral in the record: who referred and when, the board
// date, and the referral's note as written. With no referral, a single muted line in its place.
export default function ReferralSection({ referral }: { referral: TumorBoardReferral | null }) {
  if (!referral) {
    const { title, icon: Icon } = CASE_SECTIONS.referral;
    return (
      <section
        id={sectionAnchor('referral')}
        aria-label={title}
        style={{ scrollMarginTop: 'var(--case-scroll-offset, 8rem)' }}
        className="flex items-center gap-2.5 rounded-2xl border border-dashed border-slate-300 bg-white px-5 py-3"
      >
        <span className="flex h-7 w-7 shrink-0 items-center justify-center rounded-lg bg-slate-100 text-slate-500" aria-hidden="true">
          <Icon size={14} strokeWidth={2.25} />
        </span>
        <span className="text-sm text-slate-600">No tumor board referral in the record.</span>
      </section>
    );
  }

  const { requester, authoredOn, boardDate, notes } = referral;
  const referredBy = requester || authoredOn ? (
    <span>
      {requester ? <>Referred by <span className="font-semibold text-slate-900">{requester}</span></> : 'Referred'}
      {authoredOn && <> · {authoredOn}</>}
    </span>
  ) : null;

  return (
    <SectionCard
      section="referral"
      actions={(referredBy || boardDate) && (
        <div className="flex flex-wrap items-center gap-2 text-sm text-slate-600">
          {referredBy}
          {boardDate && (
            <span className="whitespace-nowrap rounded-full bg-slate-900 px-2.5 py-0.5 text-xs font-medium text-white">
              Tumor board {boardDate}
            </span>
          )}
        </div>
      )}
    >
      {notes.length > 0 ? (
        <div className="flex max-w-5xl flex-col gap-2 text-[15px] leading-relaxed text-slate-800">
          {notes.map((note, position) => <p key={position}>{note}</p>)}
        </div>
      ) : (
        <p className="text-sm text-slate-500">The referral has no note.</p>
      )}
    </SectionCard>
  );
}
