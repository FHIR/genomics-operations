import {
  AllergyIntolerance,
  CodeableConcept,
  Condition,
  DiagnosticReport,
  FamilyMemberHistory,
  FhirResource,
  Observation,
  ObservationValue,
  Reference,
} from '@/types/fhir';
import { CrosswalkRow } from '@/cancerFilter/tumorTypeCrosswalk';
import { hasTimelineDate, isMedication } from '@/lib/treatmentTimeline';
import { ClinicalRecord } from '@/services/clinicalService';
import {
  CODE_SYSTEM_NAMES,
  conceptLabel,
  DISEASE_STATUS_CODES,
  effectiveDate,
  firstCode,
  formatDate,
  getPrimaryCancerConditions,
  hasCode,
  hasProfile,
  isEnteredInError,
  isSecondaryCancer,
  latestFirst,
  LOINC,
  SNOMED,
  STAGE_GROUP_CODES,
} from '@/lib/clinicalSummary';

// Sorts a patient's clinical record into the Case summary sections (mtb-clinical-data-design.md 4.2, 4.3).
// Staging, TNM categories and disease status are recognized first, by mCODE profile or code, because
// authors file them under different Observation categories (exam, imaging, laboratory). Everything else
// goes by Observation.category and resource type, and whatever isn't recognized lands in Other resources.

const TNM_CATEGORY_PROFILES = [
  'mcode-tnm-primary-tumor-category',
  'mcode-tnm-regional-nodes-category',
  'mcode-tnm-distant-metastases-category',
];
const TNM_CATEGORY_CODES = [
  `${SNOMED}|399504009`, // cT category
  `${SNOMED}|399534004`, // cN category
  `${SNOMED}|399387003`, // cM category
  `${SNOMED}|384625004`, // pT category
  `${SNOMED}|371494008`, // pN category
  `${SNOMED}|371497001`, // pM category
  `${LOINC}|21905-5`, // Primary tumor.clinical
  `${LOINC}|21906-3`, // Regional lymph nodes.clinical
  `${LOINC}|21907-1`, // Distant metastases.clinical
  `${LOINC}|21899-0`, // Primary tumor.pathology
  `${LOINC}|21900-6`, // Regional lymph nodes.pathology
  `${LOINC}|21901-4`, // Distant metastases.pathology
];
const HISTOLOGY_EXTENSION = 'http://hl7.org/fhir/us/mcode/StructureDefinition/mcode-histology-morphology-behavior';
const LATERALITY_EXTENSIONS = [
  'http://hl7.org/fhir/us/mcode/StructureDefinition/mcode-laterality-qualifier',
  'http://hl7.org/fhir/us/mcode/StructureDefinition/mcode-laterality',
];

export interface ObservationGroup {
  key: string;
  label: string;
  // Newest first, exact duplicates removed
  observations: Observation[];
}

export interface StageGroup {
  observation: Observation;
  members: Observation[];
}

export interface PrimaryCancer {
  condition: Condition;
  stageGroups: StageGroup[];
}

export interface ReportGroup {
  label: string;
  // Newest first
  reports: DiagnosticReport[];
}

export interface ComorbidityNote {
  date?: string;
  labels: string[];
}

export interface CaseSummaryData {
  primaryCancers: PrimaryCancer[];
  // Stage groups and TNM categories not tied to a recognized primary cancer
  otherStaging: StageGroup[];
  looseTnmCategories: Observation[];
  metastases: Condition[];
  labs: ObservationGroup[];
  imaging: ObservationGroup[];
  assessments: ObservationGroup[];
  vitals: ObservationGroup[];
  otherConditions: Condition[];
  // From comorbidity observations: what they add beyond the Conditions already listed, by assessment date
  comorbiditiesAbsent: ComorbidityNote[];
  comorbiditiesNotListed: ComorbidityNote[];
  allergies: AllergyIntolerance[];
  reportGroups: ReportGroup[];
  familyHistory: FamilyMemberHistory[];
  // For the treatment timeline: dated medications and procedures, and disease status about the main cancer
  medications: FhirResource[];
  procedures: FhirResource[];
  diseaseStatuses: Observation[];
  index: Map<string, FhirResource>;
}

export function resourceKey(resource: FhirResource) {
  return `${resource.resourceType}/${resource.id}`;
}

// "Condition/abc" or an absolute URL ending in it
export function resolveReference(reference: Reference | undefined, index: Map<string, FhirResource>) {
  const parts = reference?.reference?.split('/') ?? [];
  return parts.length >= 2 ? index.get(parts.slice(-2).join('/')) : undefined;
}

function hasCategory(observation: Observation, ...codes: string[]) {
  return observation.category?.some((category) => category.coding?.some((coding) => codes.includes(coding.code ?? ''))) ?? false;
}

function codeKey(concept: CodeableConcept | undefined) {
  const coding = concept?.coding?.[0];
  return coding?.code ? `${coding.system}|${coding.code}` : `text|${concept?.text ?? ''}`;
}

function conditionDate(condition: Condition) {
  return condition.onsetDateTime ?? condition.onsetPeriod?.start ?? condition.recordedDate;
}

export function reportDate(report: DiagnosticReport) {
  return report.effectiveDateTime ?? report.effectivePeriod?.start ?? report.issued;
}

function formatQuantity(quantity: { value?: number; unit?: string; code?: string } | undefined) {
  if (quantity?.value === undefined) {
    return undefined;
  }
  const unit = quantity.unit ?? quantity.code;
  if (!unit) {
    return String(quantity.value);
  }
  return unit === '%' ? `${quantity.value}%` : `${quantity.value} ${unit}`;
}

function valueText(value: ObservationValue): string | undefined {
  if (value.valueQuantity) {
    return formatQuantity(value.valueQuantity);
  }
  if (value.valueCodeableConcept) {
    return conceptLabel(value.valueCodeableConcept);
  }
  if (value.valueString !== undefined) {
    return value.valueString;
  }
  if (value.valueInteger !== undefined) {
    return String(value.valueInteger);
  }
  if (value.valueBoolean !== undefined) {
    return value.valueBoolean ? 'Yes' : 'No';
  }
  if (value.valueRange) {
    return [formatQuantity(value.valueRange.low), formatQuantity(value.valueRange.high)].filter(Boolean).join(' – ');
  }
  if (value.valueRatio) {
    return `${formatQuantity(value.valueRatio.numerator) ?? '?'} / ${formatQuantity(value.valueRatio.denominator) ?? '?'}`;
  }
  if (value.valueDateTime) {
    return formatDate(value.valueDateTime);
  }
  if (value.dataAbsentReason) {
    return `Not available (${conceptLabel(value.dataAbsentReason)})`;
  }
  return undefined;
}

/**
 * The observation's value as text: value[x], else its components, else the conditions its extensions
 * point to (e.g. mCODE comorbidities). Undefined when it has none of these.
 */
export function observationValueText(observation: Observation, index: Map<string, FhirResource>): string | undefined {
  const own = valueText(observation);
  if (own) {
    return own;
  }
  // A component's name is shown only when it has readable text; a bare code adds nothing for the reader
  const components = (observation.component ?? [])
    .map((component) => {
      const text = valueText(component);
      const hasName = Boolean(component.code?.text || component.code?.coding?.some((coding) => coding.display));
      return text && hasName ? `${conceptLabel(component.code)}: ${text}` : text;
    })
    .filter(Boolean);
  if (components.length > 0) {
    return components.join('; ');
  }
  const referenced = (observation.extension ?? [])
    .map((extension) => {
      const target = resolveReference(extension.valueReference, index) as Condition | undefined;
      const label = conceptLabel(target?.code) ?? extension.valueReference?.display ?? conceptLabel(extension.valueCodeableConcept);
      return label && /absent/i.test(extension.url) ? `no ${label}` : label;
    })
    .filter(Boolean);
  return referenced.length > 0 ? referenced.join(', ') : undefined;
}

/**
 * The numeric value used for trend lines, with its unit; undefined for non-numeric results.
 */
export function numericValue(observation: Observation): { value: number; unit: string } | undefined {
  if (observation.valueQuantity?.value !== undefined) {
    return { value: observation.valueQuantity.value, unit: observation.valueQuantity.code ?? observation.valueQuantity.unit ?? '' };
  }
  if (observation.valueInteger !== undefined) {
    return { value: observation.valueInteger, unit: '' };
  }
  return undefined;
}

export function referenceRangeText(observation: Observation) {
  const range = observation.referenceRange?.[0];
  if (!range) {
    return undefined;
  }
  if (range.text) {
    return range.text;
  }
  const low = range.low?.value;
  const high = range.high?.value;
  if (low !== undefined && high !== undefined) {
    return `${low}–${high}`;
  }
  if (low !== undefined) {
    return `≥${low}`;
  }
  return high !== undefined ? `≤${high}` : undefined;
}

const HIGH_FLAGS = ['H', 'HH', 'HU', '>'];
const LOW_FLAGS = ['L', 'LL', 'LU', '<'];

export function interpretationFlag(observation: Observation): { text: string; direction?: 'high' | 'low' } | undefined {
  const interpretation = observation.interpretation?.[0];
  const code = interpretation?.coding?.find((coding) => coding.system?.includes('ObservationInterpretation'))?.code;
  if (code === 'N') {
    return undefined;
  }
  if (code && HIGH_FLAGS.includes(code)) {
    return { text: `▲ ${code}`, direction: 'high' };
  }
  if (code && LOW_FLAGS.includes(code)) {
    return { text: `▼ ${code}`, direction: 'low' };
  }
  const label = conceptLabel(interpretation);
  return label ? { text: label } : undefined;
}

export function bodySiteText(condition: Condition) {
  return (condition.bodySite ?? []).map((site) => {
    const laterality = site.extension
      ?.filter((extension) => LATERALITY_EXTENSIONS.includes(extension.url))
      .map((extension) => conceptLabel(extension.valueCodeableConcept))
      .filter(Boolean);
    const label = conceptLabel(site) ?? 'Unnamed site';
    return laterality?.length ? `${label} (${laterality.join(', ')})` : label;
  });
}

export function histologyText(condition: Condition) {
  const histology = condition.extension?.find((extension) => extension.url === HISTOLOGY_EXTENSION);
  return conceptLabel(histology?.valueCodeableConcept);
}

// Every coding of a concept, e.g. "SNOMED 254626006 Adenocarcinoma of lung · ICD-10-CM C34.31"
export function codingsText(concept: CodeableConcept | undefined) {
  return (concept?.coding ?? [])
    .map((coding) => [CODE_SYSTEM_NAMES[coding.system ?? ''] ?? coding.system, coding.code, coding.display].filter(Boolean).join(' '))
    .join(' · ');
}

function isStageGroup(observation: Observation) {
  return hasProfile(observation, 'mcode-tnm-stage-group') || hasCode(observation.code, STAGE_GROUP_CODES);
}

function isTnmCategory(observation: Observation) {
  return TNM_CATEGORY_PROFILES.some((profile) => hasProfile(observation, profile)) || hasCode(observation.code, TNM_CATEGORY_CODES);
}

const COMORBIDITY_CODES = [`${SNOMED}|398192003`]; // Co-morbid conditions

function isComorbidities(observation: Observation) {
  return hasProfile(observation, 'mcode-comorbidities') || hasCode(observation.code, COMORBIDITY_CODES);
}

function isDiseaseStatus(observation: Observation) {
  return hasProfile(observation, 'mcode-cancer-disease-status') || hasCode(observation.code, DISEASE_STATUS_CODES);
}

// Group observations by code, newest first; exact repeats (same code, date and value) are dropped
function groupObservations(observations: Observation[], index: Map<string, FhirResource>): ObservationGroup[] {
  const groups = new Map<string, ObservationGroup>();
  const readableKeys = new Set<string>();
  for (const observation of [...observations].sort(latestFirst)) {
    const key = codeKey(observation.code);
    const hasReadableLabel = Boolean(observation.code?.text || observation.code?.coding?.some((coding) => coding.display));
    let group = groups.get(key);
    if (!group) {
      group = { key, label: conceptLabel(observation.code) ?? 'Unnamed observation', observations: [] };
      groups.set(key, group);
      if (hasReadableLabel) {
        readableKeys.add(key);
      }
    } else if (hasReadableLabel && !readableKeys.has(key)) {
      // Prefer a human-readable name from any result over a bare code from the newest one
      group.label = conceptLabel(observation.code) ?? group.label;
      readableKeys.add(key);
    }
    const signature = `${effectiveDate(observation)}|${observationValueText(observation, index)}`;
    if (!group.observations.some((existing) => `${effectiveDate(existing)}|${observationValueText(existing, index)}` === signature)) {
      group.observations.push(observation);
    }
  }
  return [...groups.values()].sort((left, right) => left.label.localeCompare(right.label));
}

export function classifyClinicalRecord(record: ClinicalRecord, crosswalk: CrosswalkRow[]): CaseSummaryData {
  const index = new Map(record.resources.map((resource) => [resourceKey(resource), resource]));
  // Resources already placed in a section (staging is placed before observations are sorted by category)
  const placed = new Set<string>();
  const place = (resource: FhirResource) => placed.add(resourceKey(resource));

  // The Patient is shown in the banner; entered-in-error records and cancelled observations aren't shown
  const usable = record.resources.filter((resource) => resource.resourceType !== 'Patient'
    && !isEnteredInError(resource)
    && !(resource.resourceType === 'Observation' && resource.status === 'cancelled'));
  const ofType = <T extends FhirResource>(type: string) => usable.filter((resource): resource is T => resource.resourceType === type);

  // Cancer conditions
  const conditions = ofType<Condition>('Condition');
  const primaries = getPrimaryCancerConditions(conditions, crosswalk);
  const primaryKeys = new Set(primaries.map(resourceKey));
  const metastases = conditions
    .filter((condition) => !primaryKeys.has(resourceKey(condition)) && isSecondaryCancer(condition))
    .sort((left, right) => (conditionDate(right) ?? '').localeCompare(conditionDate(left) ?? ''));
  const metastasisKeys = new Set(metastases.map(resourceKey));
  const otherConditions = conditions
    .filter((condition) => !primaryKeys.has(resourceKey(condition)) && !metastasisKeys.has(resourceKey(condition)))
    .sort((left, right) => (conceptLabel(left.code) ?? '').localeCompare(conceptLabel(right.code) ?? ''));
  [...primaries, ...metastases].forEach(place);
  otherConditions.forEach(place);

  // Staging: stage groups with their T/N/M members, attached to the cancer they're about
  const observations = ofType<Observation>('Observation');
  const stageGroups: StageGroup[] = observations.filter(isStageGroup).sort(latestFirst).map((observation) => ({
    observation,
    members: (observation.hasMember ?? [])
      .map((member) => resolveReference(member, index))
      .filter((member): member is Observation => member?.resourceType === 'Observation'),
  }));
  const memberKeys = new Set(stageGroups.flatMap((group) => group.members.map(resourceKey)));
  const looseTnmCategories = observations.filter((observation) => isTnmCategory(observation) && !memberKeys.has(resourceKey(observation)));
  stageGroups.forEach((group) => {
    place(group.observation);
    group.members.forEach(place);
  });
  looseTnmCategories.forEach(place);

  const focusKeys = (observation: Observation) => (observation.focus ?? [])
    .map((focus) => resolveReference(focus, index))
    .filter((target): target is FhirResource => target !== undefined)
    .map(resourceKey);
  const primaryCancers: PrimaryCancer[] = primaries.map((condition, position) => ({
    condition,
    // A stage group without a focus belongs to the main (first) cancer
    stageGroups: stageGroups.filter((group) => {
      const keys = focusKeys(group.observation);
      return keys.includes(resourceKey(condition)) || (position === 0 && !group.observation.focus?.length);
    }),
  }));
  const attached = new Set(primaryCancers.flatMap((cancer) => cancer.stageGroups.map((group) => resourceKey(group.observation))));
  const otherStaging = stageGroups.filter((group) => !attached.has(resourceKey(group.observation)));

  // Remaining observations by category
  const byCategory: Record<'labs' | 'imaging' | 'assessments' | 'vitals', Observation[]> = { labs: [], imaging: [], assessments: [], vitals: [] };
  const comorbidityObservations: Observation[] = [];
  const diseaseStatusObservations: Observation[] = [];
  for (const observation of observations) {
    if (placed.has(resourceKey(observation))) {
      continue;
    }
    if (isDiseaseStatus(observation)) {
      place(observation);
      diseaseStatusObservations.push(observation);
    } else if (isComorbidities(observation)) {
      // Its conditions are listed under Other conditions; only what it adds is shown there
      place(observation);
      comorbidityObservations.push(observation);
    } else if (hasCategory(observation, 'laboratory')) {
      byCategory.labs.push(observation);
    } else if (hasCategory(observation, 'imaging')) {
      byCategory.imaging.push(observation);
    } else if (hasCategory(observation, 'survey', 'exam')) {
      byCategory.assessments.push(observation);
    } else if (hasCategory(observation, 'vital-signs', 'social-history')) {
      byCategory.vitals.push(observation);
    }
  }
  const grouped = (section: keyof typeof byCategory) => {
    byCategory[section].forEach(place);
    return groupObservations(byCategory[section], index);
  };
  // Comorbidities: conditions recorded as absent, and present ones that aren't Conditions in this record
  const listedConditionKeys = new Set(conditions.map(resourceKey));
  const comorbiditiesAbsent: ComorbidityNote[] = [];
  const comorbiditiesNotListed: ComorbidityNote[] = [];
  for (const observation of [...comorbidityObservations].sort(latestFirst)) {
    const absent = new Set<string>();
    const notListed = new Set<string>();
    for (const extension of observation.extension ?? []) {
      const target = resolveReference(extension.valueReference, index);
      if (target && listedConditionKeys.has(resourceKey(target))) {
        continue;
      }
      const label = conceptLabel((target as Condition | undefined)?.code) ?? extension.valueReference?.display ?? conceptLabel(extension.valueCodeableConcept);
      if (label) {
        (/absent/i.test(extension.url) ? absent : notListed).add(label);
      }
    }
    const date = formatDate(effectiveDate(observation));
    if (absent.size > 0) {
      comorbiditiesAbsent.push({ date, labels: [...absent].sort() });
    }
    if (notListed.size > 0) {
      comorbiditiesNotListed.push({ date, labels: [...notListed].sort() });
    }
  }

  const labs = grouped('labs');
  const imaging = grouped('imaging');
  const assessments = grouped('assessments');
  const vitals = grouped('vitals');

  // Treatment timeline. Disease status follows the banner chip: final results about the main cancer (or with no focus)
  const medications = usable.filter((resource) => isMedication(resource) && hasTimelineDate(resource));
  const procedures = usable.filter((resource) => resource.resourceType === 'Procedure' && hasTimelineDate(resource));
  [...medications, ...procedures].forEach(place);
  const mainCancer = primaries[0];
  const diseaseStatuses = diseaseStatusObservations.filter((observation) => ['final', 'amended', 'corrected'].includes(observation.status ?? '')
    && (!mainCancer || !observation.focus?.length || focusKeys(observation).includes(resourceKey(mainCancer))));

  // Allergies: active first, then the rest; each block alphabetical
  const allergyRank = (allergy: AllergyIntolerance) => (['inactive', 'resolved'].includes(firstCode(allergy.clinicalStatus) ?? '') ? 1 : 0);
  const allergies = ofType<AllergyIntolerance>('AllergyIntolerance').sort((left, right) =>
    allergyRank(left) - allergyRank(right) || (conceptLabel(left.code) ?? '').localeCompare(conceptLabel(right.code) ?? ''));
  allergies.forEach(place);

  // Reports grouped by category label, newest first
  const reportsByLabel = new Map<string, DiagnosticReport[]>();
  for (const report of ofType<DiagnosticReport>('DiagnosticReport')) {
    const label = report.category?.map((category) => conceptLabel(category)).find(Boolean) ?? 'Uncategorized';
    reportsByLabel.set(label, [...(reportsByLabel.get(label) ?? []), report]);
    place(report);
  }
  const reportGroups = [...reportsByLabel.entries()]
    .map(([label, reports]) => ({ label, reports: reports.sort((left, right) => (reportDate(right) ?? '').localeCompare(reportDate(left) ?? '')) }))
    .sort((left, right) => left.label.localeCompare(right.label));

  const familyHistory = ofType<FamilyMemberHistory>('FamilyMemberHistory');

  return {
    primaryCancers,
    otherStaging,
    looseTnmCategories,
    metastases,
    labs,
    imaging,
    assessments,
    vitals,
    otherConditions,
    comorbiditiesAbsent,
    comorbiditiesNotListed,
    allergies,
    reportGroups,
    medications,
    procedures,
    diseaseStatuses,
    familyHistory,
    index,
  };
}
