import {
  CalendarRange,
  ClipboardList,
  FileText,
  FlaskConical,
  HeartPulse,
  LucideIcon,
  Ribbon,
  ScanLine,
  Stethoscope,
  TriangleAlert,
  Users,
} from 'lucide-react';

// The Case summary's sections, in page order: one place for each card's title, icon and color, shared by
// the cards and the jump bar so they never drift apart

// Each section's icon tile has its own color: a soft tint behind a stronger icon
export type IconAccent = 'violet' | 'blue' | 'teal' | 'indigo' | 'amber' | 'rose' | 'sky' | 'red' | 'emerald' | 'slate';

export const ACCENT_CLASSES: Record<IconAccent, string> = {
  violet: 'bg-violet-100 text-violet-700',
  blue: 'bg-blue-100 text-blue-700',
  teal: 'bg-teal-100 text-teal-700',
  indigo: 'bg-indigo-100 text-indigo-700',
  amber: 'bg-amber-100 text-amber-700',
  rose: 'bg-rose-100 text-rose-700',
  sky: 'bg-sky-100 text-sky-700',
  red: 'bg-red-100 text-red-700',
  emerald: 'bg-emerald-100 text-emerald-700',
  slate: 'bg-slate-200 text-slate-700',
};

export interface CaseSection {
  title: string;
  // Label in the jump bar
  short: string;
  icon: LucideIcon;
  accent: IconAccent;
}

export const CASE_SECTIONS = {
  diagnosis: { title: 'Cancer diagnosis & staging', short: 'Diagnosis', icon: Ribbon, accent: 'violet' },
  timeline: { title: 'Treatment timeline', short: 'Timeline', icon: CalendarRange, accent: 'blue' },
  labs: { title: 'Lab & pathology results', short: 'Labs', icon: FlaskConical, accent: 'teal' },
  imaging: { title: 'Imaging observations', short: 'Imaging', icon: ScanLine, accent: 'indigo' },
  assessments: { title: 'Assessments', short: 'Assessments', icon: ClipboardList, accent: 'amber' },
  vitals: { title: 'Vitals & social history', short: 'Vitals', icon: HeartPulse, accent: 'rose' },
  conditions: { title: 'Other conditions', short: 'Conditions', icon: Stethoscope, accent: 'sky' },
  allergies: { title: 'Allergies', short: 'Allergies', icon: TriangleAlert, accent: 'red' },
  family: { title: 'Family history', short: 'Family', icon: Users, accent: 'emerald' },
  reports: { title: 'Reports', short: 'Reports', icon: FileText, accent: 'slate' },
} satisfies Record<string, CaseSection>;

export type CaseSectionKey = keyof typeof CASE_SECTIONS;

export const CASE_SECTION_ORDER = Object.keys(CASE_SECTIONS) as CaseSectionKey[];

// The element id a section's card carries, for the jump bar to scroll to
export function sectionAnchor(key: CaseSectionKey) {
  return `case-section-${key}`;
}
