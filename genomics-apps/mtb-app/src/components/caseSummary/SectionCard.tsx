import { ReactNode } from 'react';

interface SectionCardProps {
  title: string;
  actions?: ReactNode;
  children: ReactNode;
  className?: string;
}

export function SectionCard({ title, actions, children, className = '' }: SectionCardProps) {
  return (
    <section aria-label={title} className={`rounded-xl border border-slate-300 bg-white p-5 shadow-sm ${className}`}>
      <div className="mb-3 flex flex-wrap items-start justify-between gap-3">
        <h2 className="text-lg font-semibold text-slate-900">{title}</h2>
        {actions}
      </div>
      {children}
    </section>
  );
}
