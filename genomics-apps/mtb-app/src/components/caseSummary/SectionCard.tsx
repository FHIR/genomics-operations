import { ReactNode } from 'react';
import { ACCENT_CLASSES, CASE_SECTIONS, CaseSectionKey, sectionAnchor } from './sections';

interface SectionCardProps {
  // Which Case summary section this is: sets the title, icon, color and jump-bar anchor
  section: CaseSectionKey;
  actions?: ReactNode;
  children: ReactNode;
  className?: string;
  // Cap the body's height (about 10 lines) and scroll inside it, so large records don't stretch the page
  scrollBody?: boolean;
}

export function SectionCard({ section, actions, children, className = '', scrollBody = false }: SectionCardProps) {
  const { title, icon: Icon, accent } = CASE_SECTIONS[section];
  return (
    <section
      id={sectionAnchor(section)}
      aria-label={title}
      // Leaves room for the sticky banner and jump bar when scrolled to
      style={{ scrollMarginTop: 'var(--case-scroll-offset, 8rem)' }}
      className={`rounded-2xl border border-slate-200 bg-white p-5 shadow-[0_1px_3px_rgba(15,23,42,0.06),0_4px_12px_rgba(15,23,42,0.04)] ${className}`}>
      <div className="mb-4 flex flex-wrap items-center justify-between gap-3 border-b border-slate-100 pb-3">
        <div className="flex items-center gap-2.5">
          <span className={`flex h-8 w-8 shrink-0 items-center justify-center rounded-lg ${ACCENT_CLASSES[accent]}`} aria-hidden="true">
            <Icon size={16} strokeWidth={2.25} />
          </span>
          <h2 className="text-base font-semibold tracking-tight text-slate-900">{title}</h2>
        </div>
        {actions}
      </div>
      {scrollBody ? <div className="max-h-64 overflow-y-auto pr-1">{children}</div> : children}
    </section>
  );
}

export type ChipTone = 'neutral' | 'dark' | 'high' | 'low' | 'danger' | 'muted';

const CHIP_TONES: Record<ChipTone, string> = {
  neutral: 'bg-slate-100 text-slate-700 ring-1 ring-inset ring-slate-200',
  dark: 'bg-slate-900 text-white',
  high: 'bg-orange-50 text-orange-800 ring-1 ring-inset ring-orange-200',
  low: 'bg-blue-50 text-blue-800 ring-1 ring-inset ring-blue-200',
  danger: 'bg-red-50 text-red-800 ring-1 ring-inset ring-red-200',
  muted: 'bg-white text-slate-500 ring-1 ring-inset ring-slate-200',
};

// A small rounded label, e.g. a stage, a body site or a lab flag
export function Chip({ children, tone = 'neutral', title }: { children: ReactNode; tone?: ChipTone; title?: string }) {
  return (
    <span title={title} className={`inline-flex items-center whitespace-nowrap rounded-full px-2.5 py-0.5 text-xs font-medium ${CHIP_TONES[tone]}`}>
      {children}
    </span>
  );
}
