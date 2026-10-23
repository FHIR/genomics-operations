import { CodeableConcept, Condition, FhirResource, Observation, Period } from '@/types/fhir';
import {
  conceptLabel,
  conditionOnset,
  effectiveDate,
  formatDate,
  hasCode,
  hasProfile,
  SNOMED,
} from '@/lib/clinicalSummary';
import { bodySiteText, CaseSummaryData, resolveReference } from '@/lib/caseSummary';

// Builds the Case summary's treatment timeline (mtb-clinical-data-design.md 4.2, "Treatment timeline rules").
// Every mark comes from a FHIR resource; nothing is inferred. Row order is fixed: active medications A–Z,
// inactive medications A–Z, radiotherapy, procedures, disease status, diagnosis & metastases.

const MCODE = 'http://hl7.org/fhir/us/mcode/StructureDefinition/';
const RADIOTHERAPY_CODES = [
  `${SNOMED}|1217123003`, // Radiotherapy course of treatment
  `${SNOMED}|108290001`, // Radiation oncology AND/OR radiotherapy
  `${SNOMED}|53438000`, // Radiation therapy procedure or service
];
const MEDICATION_TYPES = ['MedicationRequest', 'MedicationAdministration', 'MedicationStatement'];
const ACTIVE_MEDICATION_STATUSES = ['active', 'in-progress', 'intended', 'on-hold'];
// Events this close together share one marker label
const GROUPING_DAYS = 30;

// mCODE's own disease status terms (never RECIST labels)
const DISEASE_STATUS_TERMS: Record<string, { label: string; tone: MarkerTone }> = {
  '268910001': { label: '▲ Improved', tone: 'improving' },
  '359746009': { label: '● Stable', tone: 'stable' },
  '271299001': { label: '▼ Worsened', tone: 'progression' },
  '709137006': { label: 'Undetermined', tone: 'other' },
};

export type LaneGroup = 'Active medications' | 'Inactive medications' | 'Radiotherapy & procedures' | 'Disease status' | 'Diagnosis & metastases';
export type MarkerTone = 'improving' | 'stable' | 'progression' | 'other' | 'procedure' | 'diagnosis';

export interface TimelineBar {
  start: string;
  end: string;
  kind: 'medication' | 'radiotherapy';
  label: string;
  tooltip: string;
}

export interface TimelineMarker {
  date: string;
  symbol: string;
  label?: string;
  tooltip: string;
  tone: MarkerTone;
  emphasized?: boolean;
}

export interface TimelineLane {
  key: string;
  group: LaneGroup;
  // Empty for a lane that is the only one in its group; the group heading names it
  label: string;
  bars: TimelineBar[];
  markers: TimelineMarker[];
}

export interface TreatmentTimeline {
  start: string;
  end: string;
  lanes: TimelineLane[];
}

function day(value: string) {
  return value.slice(0, 10);
}

function daysBetween(left: string, right: string) {
  return Math.abs(Date.parse(day(right)) - Date.parse(day(left))) / 864e5;
}

function monthLabel(value: string) {
  return new Date(`${day(value)}T00:00:00`).toLocaleDateString('en-US', { month: 'short', year: 'numeric' });
}

export function isRadiotherapy(procedure: FhirResource) {
  const code = procedure.code as CodeableConcept | undefined;
  return hasProfile(procedure, 'mcode-radiotherapy-course-summary') || hasCode(code, RADIOTHERAPY_CODES);
}

function procedureDates(procedure: FhirResource): { start?: string; end?: string } {
  const period = procedure.performedPeriod as Period | undefined;
  const dateTime = procedure.performedDateTime as string | undefined;
  return { start: dateTime ?? period?.start, end: dateTime ?? period?.end };
}

// The medication's dates: dosing period (request), effective period or date (administration, statement)
function medicationDates(resource: FhirResource): { start?: string; end?: string; hasPeriod: boolean } {
  if (resource.resourceType === 'MedicationRequest') {
    const dosage = resource.dosageInstruction as { timing?: { repeat?: { boundsPeriod?: Period } } }[] | undefined;
    const bounds = dosage?.map((instruction) => instruction.timing?.repeat?.boundsPeriod).find((period) => period?.start);
    return bounds ? { start: bounds.start, end: bounds.end, hasPeriod: true } : { start: resource.authoredOn as string | undefined, hasPeriod: false };
  }
  const period = resource.effectivePeriod as Period | undefined;
  const dateTime = resource.effectiveDateTime as string | undefined;
  return { start: dateTime ?? period?.start, end: dateTime ?? period?.end, hasPeriod: Boolean(dateTime ?? period?.start) };
}

/**
 * Whether a medication or procedure can be placed on the timeline (it has a date).
 */
export function hasTimelineDate(resource: FhirResource) {
  if (MEDICATION_TYPES.includes(resource.resourceType)) {
    return Boolean(medicationDates(resource).start);
  }
  return resource.resourceType === 'Procedure' && Boolean(procedureDates(resource).start);
}

export function isMedication(resource: FhirResource) {
  return MEDICATION_TYPES.includes(resource.resourceType);
}

function medicationConcept(resource: FhirResource, index: Map<string, FhirResource>) {
  const concept = resource.medicationCodeableConcept as CodeableConcept | undefined;
  if (concept) {
    return concept;
  }
  const medication = resolveReference(resource.medicationReference as { reference?: string } | undefined, index);
  return medication?.code as CodeableConcept | undefined;
}

function medicationLanes(medications: FhirResource[], index: Map<string, FhirResource>, today: string): TimelineLane[] {
  const lanes = new Map<string, { label: string; isActive: boolean; bars: TimelineBar[] }>();
  for (const resource of medications) {
    const concept = medicationConcept(resource, index);
    const label = conceptLabel(concept) ?? 'Medication (no code)';
    // Records for the same drug share a lane
    const coding = concept?.coding?.[0];
    const key = coding?.code ? `${coding.system}|${coding.code}` : `text|${label}`;
    const status = resource.status as string | undefined;
    const isActive = ACTIVE_MEDICATION_STATUSES.includes(status ?? '');
    const statusReason = conceptLabel(resource.statusReason as CodeableConcept | undefined);
    const dosageText = (resource.dosageInstruction as { text?: string }[] | undefined)?.map((instruction) => instruction.text).filter(Boolean).join('; ');
    const { start, end, hasPeriod } = medicationDates(resource);
    if (!start) {
      continue;
    }

    const statusLabel = [status, statusReason].filter(Boolean).join(' · ');
    const bar: TimelineBar = hasPeriod
      ? {
        start,
        end: end ?? (isActive ? today : start),
        kind: 'medication',
        label: `${formatDate(start)} → ${end ? formatDate(end) : isActive ? 'ongoing' : '?'}${statusLabel ? ` · ${statusLabel}` : ''}`,
        tooltip: [label, `${formatDate(start)} → ${end ? formatDate(end) : isActive ? 'ongoing' : 'no end date'}`, statusLabel, dosageText].filter(Boolean).join('\n'),
      }
      : {
        start,
        end: isActive ? today : start,
        kind: 'medication',
        // Drawn like any other bar; the label says the dosing period isn't recorded
        label: `ordered ${formatDate(start)}, no dosing period${statusLabel ? ` · ${statusLabel}` : ''}`,
        tooltip: [label, `Ordered ${formatDate(start)}; no dosing period recorded`, statusLabel, dosageText].filter(Boolean).join('\n'),
      };

    const lane = lanes.get(key) ?? { label, isActive: false, bars: [] };
    lane.isActive = lane.isActive || isActive;
    lane.bars.push(bar);
    lanes.set(key, lane);
  }

  const toLane = ([key, lane]: [string, { label: string; isActive: boolean; bars: TimelineBar[] }]): TimelineLane => ({
    key: `medication:${key}`,
    group: lane.isActive ? 'Active medications' : 'Inactive medications',
    label: lane.label,
    bars: lane.bars.sort((left, right) => left.start.localeCompare(right.start)),
    markers: [],
  });
  // Alphabetical by drug name: a leading volume such as "10 ML" is skipped for sorting (not for display)
  const sortKey = (label: string) => label.replace(/^[\d.]+\s*(ML|MG|MCG|G|L)\s+/i, '');
  const byLabel = (left: TimelineLane, right: TimelineLane) => sortKey(left.label).localeCompare(sortKey(right.label), undefined, { sensitivity: 'base' });
  const all = [...lanes.entries()].map(toLane);
  return [
    ...all.filter((lane) => lane.group === 'Active medications').sort(byLabel),
    ...all.filter((lane) => lane.group === 'Inactive medications').sort(byLabel),
  ];
}

// Body site, dose and fractions (mCODE dose extension, when present) and status
function radiotherapyLabel(procedure: FhirResource) {
  const sites = bodySiteText(procedure as Condition);
  const doses = (procedure.extension ?? [])
    .filter((extension) => extension.url === `${MCODE}mcode-radiotherapy-dose-delivered-to-volume`)
    .map((extension) => {
      const parts = (extension.extension ?? []) as { url: string; valueQuantity?: { value?: number; unit?: string }; valueUnsignedInt?: number }[];
      const dose = parts.find((part) => part.url === 'totalDoseDelivered')?.valueQuantity as { value?: number; unit?: string; code?: string } | undefined;
      const fractions = parts.find((part) => part.url === 'fractionsDelivered')?.valueUnsignedInt;
      return [dose?.value !== undefined ? `${dose.value} ${dose.unit ?? dose.code ?? ''}`.trim() : undefined, fractions !== undefined ? `${fractions} fractions` : undefined]
        .filter(Boolean).join(' / ');
    })
    .filter(Boolean);
  return [sites.join(', '), doses.join('; '), procedure.status as string | undefined].filter(Boolean).join(' · ');
}

function radiotherapyLanes(procedures: FhirResource[], today: string): TimelineLane[] {
  return procedures
    .filter(isRadiotherapy)
    .map((procedure) => ({ procedure, ...procedureDates(procedure) }))
    .filter((course): course is { procedure: FhirResource; start: string; end?: string } => Boolean(course.start))
    .sort((left, right) => left.start.localeCompare(right.start))
    .map(({ procedure, start, end }) => {
      const label = radiotherapyLabel(procedure);
      const ongoing = !end && ['in-progress', 'preparation'].includes(procedure.status as string);
      return {
        key: `radiotherapy:${procedure.id}`,
        group: 'Radiotherapy & procedures' as const,
        label: conceptLabel(procedure.code as CodeableConcept | undefined) ?? 'Radiotherapy',
        bars: [{
          start,
          end: end ?? (ongoing ? today : start),
          kind: 'radiotherapy' as const,
          label,
          tooltip: [`${formatDate(start)} → ${end ? formatDate(end) : ongoing ? 'ongoing' : 'no end date'}`, label].join('\n'),
        }],
        markers: [],
      };
    });
}

// Group dated events that fall within GROUPING_DAYS of the first event in the group
function groupByProximity<T extends { date: string }>(events: T[]): T[][] {
  const groups: T[][] = [];
  for (const event of [...events].sort((left, right) => left.date.localeCompare(right.date))) {
    const current = groups[groups.length - 1];
    if (current && daysBetween(current[0].date, event.date) <= GROUPING_DAYS) {
      current.push(event);
    } else {
      groups.push([event]);
    }
  }
  return groups;
}

// "(2024-05-20)" for one day, "(May 2024)" within a month, "(Apr – May 2024)" across months
function groupLabel(labels: string[], dates: string[]) {
  const sorted = [...dates].sort();
  const first = sorted[0];
  const last = sorted[sorted.length - 1];
  let when: string;
  if (day(first) === day(last)) {
    when = formatDate(first) as string;
  } else if (first.slice(0, 7) === last.slice(0, 7)) {
    when = monthLabel(first);
  } else if (first.slice(0, 4) === last.slice(0, 4)) {
    when = `${monthLabel(first).split(' ')[0]} – ${monthLabel(last)}`;
  } else {
    when = `${monthLabel(first)} – ${monthLabel(last)}`;
  }
  return `${labels.join(' · ')} (${when})`;
}

function procedureLane(procedures: FhirResource[]): TimelineLane | null {
  const events = procedures
    .filter((procedure) => !isRadiotherapy(procedure))
    .map((procedure) => ({ date: procedureDates(procedure).start as string, label: conceptLabel(procedure.code as CodeableConcept | undefined) ?? 'Procedure' }))
    .filter((event) => event.date);
  if (events.length === 0) {
    return null;
  }
  const markers = groupByProximity(events).flatMap((group) => group.map((event, position) => ({
    date: event.date,
    symbol: '◆',
    tone: 'procedure' as const,
    // One label per group, on its first marker
    label: position === 0 ? groupLabel(group.map((item) => item.label), group.map((item) => item.date)) : undefined,
    tooltip: `${formatDate(event.date)} ${event.label}`,
  })));
  return { key: 'procedures', group: 'Radiotherapy & procedures', label: 'Procedures', bars: [], markers };
}

function diseaseStatusLane(statuses: Observation[]): TimelineLane | null {
  const sorted = [...statuses]
    .filter((observation) => effectiveDate(observation))
    .sort((left, right) => (effectiveDate(left) as string).localeCompare(effectiveDate(right) as string));
  if (sorted.length === 0) {
    return null;
  }
  const markers = sorted.map((observation, position) => {
    const code = observation.valueCodeableConcept?.coding?.find((coding) => coding.system === SNOMED)?.code;
    const term = code ? DISEASE_STATUS_TERMS[code] : undefined;
    const label = term?.label ?? conceptLabel(observation.valueCodeableConcept) ?? 'Status';
    const date = effectiveDate(observation) as string;
    const evidence = (observation.extension ?? [])
      .filter((extension) => /evidence-type$/.test(extension.url))
      .map((extension) => conceptLabel(extension.valueCodeableConcept))
      .filter(Boolean);
    const isLatest = position === sorted.length - 1;
    return {
      date,
      symbol: '',
      label: isLatest ? `${label} ${formatDate(date)}` : label,
      tone: term?.tone ?? 'other',
      emphasized: isLatest,
      tooltip: [`${formatDate(date)} ${conceptLabel(observation.valueCodeableConcept) ?? label}`, evidence.length ? `Evidence: ${evidence.join(', ')}` : undefined]
        .filter(Boolean).join('\n'),
    };
  });
  return { key: 'disease-status', group: 'Disease status', label: '', bars: [], markers };
}

function diagnosisLane(primaries: Condition[], metastases: Condition[]): TimelineLane | null {
  const events = [
    ...primaries.map((condition) => ({ date: conditionOnset(condition), label: 'Dx', detail: conceptLabel(condition.code) ?? 'Primary cancer' })),
    ...metastases.map((condition) => {
      const detail = conceptLabel(condition.code) ?? 'Metastasis';
      return { date: conditionOnset(condition), label: detail, detail };
    }),
  ].filter((event): event is { date: string; label: string; detail: string } => Boolean(event.date));
  if (events.length === 0) {
    return null;
  }
  const markers = groupByProximity(events).flatMap((group) => {
    const metastasisCount = group.filter((event) => event.label !== 'Dx').length;
    const hasDx = group.some((event) => event.label === 'Dx');
    // Spell out a single event; summarize a cluster ("Dx + 2 metastases")
    const labels = group.length === 1
      ? [group[0].label]
      : [hasDx ? 'Dx' : undefined, metastasisCount ? `${hasDx ? '+ ' : ''}${metastasisCount} metastas${metastasisCount === 1 ? 'is' : 'es'}` : undefined].filter((part): part is string => Boolean(part));
    return group.map((event, position) => ({
      date: event.date,
      symbol: '◆',
      tone: 'diagnosis' as const,
      label: position === 0 ? groupLabel([labels.join(' ')], group.map((item) => item.date)) : undefined,
      tooltip: group.map((item) => `${formatDate(item.date)} ${item.detail}`).join('\n'),
    }));
  });
  return { key: 'diagnosis', group: 'Diagnosis & metastases', label: '', bars: [], markers };
}

export function buildTreatmentTimeline(data: CaseSummaryData, today: string): TreatmentTimeline | null {
  const lanes = [
    ...medicationLanes(data.medications, data.index, today),
    ...radiotherapyLanes(data.procedures, today),
    procedureLane(data.procedures),
    diseaseStatusLane(data.diseaseStatuses),
    diagnosisLane(data.primaryCancers.map((cancer) => cancer.condition), data.metastases),
  ].filter((lane): lane is TimelineLane => lane !== null);
  if (lanes.length === 0) {
    return null;
  }
  const dates = lanes.flatMap((lane) => [...lane.bars.flatMap((bar) => [bar.start, bar.end]), ...lane.markers.map((marker) => marker.date)]);
  const start = dates.map(day).sort()[0];
  return { start, end: day(today) > start ? day(today) : start, lanes };
}
