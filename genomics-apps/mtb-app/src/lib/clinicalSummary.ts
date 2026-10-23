import { CodeableConcept, Condition, FhirResource, Observation, Patient } from '@/types/fhir';
import { CaseDx, CrosswalkRow, mapCodingsToTumorType } from '@/cancerFilter/tumorTypeCrosswalk';

// Selection rules for the patient banner (mtb-clinical-data-design.md sections 4.2 and 4.3). Every item
// is recognized by mCODE profile OR standard code, so data without declared profiles still shows up.

export const MCODE = 'http://hl7.org/fhir/us/mcode/StructureDefinition/';
export const SNOMED = 'http://snomed.info/sct';
export const LOINC = 'http://loinc.org';
export const ICD10CM = 'http://hl7.org/fhir/sid/icd-10-cm';

export const STAGE_GROUP_CODES = [
  `${SNOMED}|399537006`, // Clinical TNM stage grouping
  `${SNOMED}|399588009`, // Pathologic TNM stage grouping
  `${SNOMED}|399390009`, // TNM stage grouping
  `${LOINC}|21908-9`, // Stage group.clinical Cancer
  `${LOINC}|21902-2`, // Stage group.pathology Cancer
  `${LOINC}|21914-7`, // Stage group.other Cancer
];
const ECOG_CODES = [`${LOINC}|89247-1`];
export const DISEASE_STATUS_CODES = [`${LOINC}|97509-4`];
const FINAL_STATUSES = ['final', 'amended', 'corrected'];
// STU4 name, and the name used by earlier mCODE versions
const EVIDENCE_TYPE_EXTENSIONS = [`${MCODE}mcode-cancer-disease-status-evidence-type`, `${MCODE}mcode-evidence-type`];

export type DiseaseStatusTrend = 'improving' | 'stable' | 'progression' | 'other';

export interface DiseaseStatus {
  label: string;
  trend: DiseaseStatusTrend;
  date?: string;
  tooltip: string;
}

const DISEASE_STATUS_VALUES: Record<string, { label: string; trend: DiseaseStatusTrend }> = {
  '268910001': { label: '▲ Improving', trend: 'improving' },
  '359746009': { label: '● Stable', trend: 'stable' },
  '271299001': { label: '▼ Progression', trend: 'progression' },
  '709137006': { label: 'Status undetermined', trend: 'other' },
};

// Condition.clinicalStatus values worth a chip when there is no CancerDiseaseStatus; plain `active` isn't
const INFORMATIVE_CLINICAL_STATUSES: Record<string, string> = {
  recurrence: 'Recurrence',
  relapse: 'Relapse',
  remission: 'Remission',
  resolved: 'Resolved',
};

export const CODE_SYSTEM_NAMES: Record<string, string> = {
  [SNOMED]: 'SNOMED',
  [ICD10CM]: 'ICD-10-CM',
  [LOINC]: 'LOINC',
  'http://terminology.hl7.org/CodeSystem/icd-o-3': 'ICD-O-3',
  'http://www.nlm.nih.gov/research/umls/rxnorm': 'RxNorm',
};

/**
 * Label for a coded value: the author's text, else the first coding's display, else the code (with its
 * system's name when known, e.g. "LOINC 89247-1"). SNOMED displays drop their semantic tag, e.g.
 * "Liver cell carcinoma (disorder)" becomes "Liver cell carcinoma".
 */
export function conceptLabel(concept: CodeableConcept | undefined): string | undefined {
  if (concept?.text) {
    return concept.text;
  }
  const coding = concept?.coding?.find((candidate) => candidate.display) ?? concept?.coding?.[0];
  if (coding?.display) {
    return coding.system === SNOMED ? coding.display.replace(/\s+\([a-z ]+\)$/, '') : coding.display;
  }
  if (!coding?.code) {
    return undefined;
  }
  const systemName = CODE_SYSTEM_NAMES[coding.system ?? ''];
  return systemName ? `${systemName} ${coding.code}` : coding.code;
}

export function hasProfile(resource: FhirResource, profileName: string) {
  return resource.meta?.profile?.some((profile) => profile.split('|')[0] === `${MCODE}${profileName}`) ?? false;
}

export function hasCode(concept: CodeableConcept | undefined, codes: string[]) {
  return concept?.coding?.some((coding) => codes.includes(`${coding.system}|${coding.code}`)) ?? false;
}

export function firstCode(concept: CodeableConcept | undefined) {
  return concept?.coding?.[0]?.code;
}

// Records marked as entered in error, by verificationStatus (Condition, AllergyIntolerance) or status
export function isEnteredInError(resource: FhirResource) {
  const verificationStatus = resource.verificationStatus as CodeableConcept | undefined;
  return firstCode(verificationStatus) === 'entered-in-error' || resource.status === 'entered-in-error';
}

export function effectiveDate(observation: Observation) {
  return observation.effectiveDateTime ?? observation.effectivePeriod?.start;
}

// Latest first by effective date, breaking ties by issued and then lastUpdated
export function latestFirst(left: Observation, right: Observation) {
  const keys = (observation: Observation) => [effectiveDate(observation), observation.issued, observation.meta?.lastUpdated];
  const leftKeys = keys(left);
  const rightKeys = keys(right);
  for (let index = 0; index < leftKeys.length; index++) {
    const comparison = (rightKeys[index] ?? '').localeCompare(leftKeys[index] ?? '');
    if (comparison !== 0) {
      return comparison;
    }
  }
  return 0;
}

function finalObservations(resources: FhirResource[], profileName: string, codes: string[]) {
  return resources
    .filter((resource): resource is Observation => resource.resourceType === 'Observation')
    .filter((observation) => FINAL_STATUSES.includes(observation.status ?? ''))
    .filter((observation) => hasProfile(observation, profileName) || hasCode(observation.code, codes))
    .sort(latestFirst);
}

// An observation about a cancer belongs to it when its focus points at it; no focus is accepted
function isAboutCondition(observation: Observation, condition: Condition | undefined) {
  if (!condition || !observation.focus?.length) {
    return true;
  }
  return observation.focus.some((focus) => focus.reference?.endsWith(`Condition/${condition.id}`));
}

export function formatDate(value: string | undefined) {
  return value?.slice(0, 10);
}

export interface PatientSummary {
  name?: string;
  age?: number;
  sex?: string;
  birthDate?: string;
}

export function getPatientSummary(patient: Patient, today = new Date()): PatientSummary {
  const name = patient.name?.find((candidate) => candidate.use === 'official') ?? patient.name?.[0];
  const nameText = name?.text ?? [...(name?.given ?? []), name?.family].filter(Boolean).join(' ');

  let age: number | undefined;
  if (patient.birthDate && /^\d{4}-\d{2}-\d{2}/.test(patient.birthDate)) {
    const [year, month, day] = patient.birthDate.split('-').map(Number);
    age = today.getFullYear() - year;
    if (today.getMonth() + 1 < month || (today.getMonth() + 1 === month && today.getDate() < day)) {
      age -= 1;
    }
  }

  return { name: nameText || undefined, age, sex: patient.gender, birthDate: patient.birthDate };
}

// Malignant neoplasm ICD-10-CM codes, excluding secondary sites (C77-C79)
function isPrimaryMalignancyIcd10(concept: CodeableConcept | undefined) {
  return concept?.coding?.some((coding) => coding.system === ICD10CM && /^C(0\d|[1-6]\d|7[0-6]|8\d|9[0-6])/.test(coding.code ?? '')) ?? false;
}

/**
 * Primary cancer conditions, most relevant first: still active (or recurrent) before the rest, then newest.
 * Recognized by the mCODE profile, a crosswalk code or a primary-malignancy ICD-10-CM code.
 */
export function getPrimaryCancerConditions(resources: FhirResource[], crosswalk: CrosswalkRow[]): Condition[] {
  const activeStatuses = ['active', 'recurrence', 'relapse'];
  const conditionDate = (condition: Condition) => condition.onsetDateTime ?? condition.recordedDate ?? '';

  return resources
    .filter((resource): resource is Condition => resource.resourceType === 'Condition')
    .filter((condition) => !isEnteredInError(condition) && !hasProfile(condition, 'mcode-secondary-cancer-condition'))
    .filter((condition) => hasProfile(condition, 'mcode-primary-cancer-condition')
      || mapCodingsToTumorType(condition.code?.coding, crosswalk).geneListCategory !== undefined
      || isPrimaryMalignancyIcd10(condition.code))
    .sort((left, right) => {
      const leftActive = activeStatuses.includes(firstCode(left.clinicalStatus) ?? '') ? 0 : 1;
      const rightActive = activeStatuses.includes(firstCode(right.clinicalStatus) ?? '') ? 0 : 1;
      return leftActive - rightActive || conditionDate(right).localeCompare(conditionDate(left));
    });
}

/**
 * Dx dropdown entries for the patient's primary cancers, labeled as coded and mapped through the crosswalk.
 */
export function getPatientDxOptions(conditions: Condition[], crosswalk: CrosswalkRow[]): CaseDx[] {
  return conditions.map((condition) => ({
    key: `condition:${condition.id}`,
    label: conceptLabel(condition.code) ?? 'Primary cancer (no code)',
    source: 'patient',
    conditionId: condition.id,
    ...mapCodingsToTumorType(condition.code?.coding, crosswalk),
  }));
}

// Shown as a caption and a value, e.g. "Stage" "IVB (AJCC)", plus an optional muted note such as a year
export interface BannerItem {
  caption: string;
  value: string;
  note?: string;
  tooltip?: string;
}

export function isSecondaryCancer(condition: Condition) {
  return hasProfile(condition, 'mcode-secondary-cancer-condition')
    || (condition.code?.coding?.some((coding) => coding.system === ICD10CM && /^C7[7-9]/.test(coding.code ?? '')) ?? false);
}

export function conditionOnset(condition: Condition) {
  return condition.onsetDateTime ?? condition.onsetPeriod?.start ?? condition.recordedDate;
}

/**
 * Current metastases (secondary cancer conditions) with an onset after the given staging date, earliest
 * first. AJCC stage isn't reassigned when disease spreads, so these are what an older stage doesn't reflect.
 */
export function getMetastasesAfter(resources: FhirResource[], stageDate: string | undefined): Condition[] {
  if (!stageDate) {
    return [];
  }
  const day = stageDate.slice(0, 10);
  return resources
    .filter((resource): resource is Condition => resource.resourceType === 'Condition')
    .filter((condition) => isSecondaryCancer(condition) && !isEnteredInError(condition))
    .filter((condition) => !['inactive', 'resolved', 'remission'].includes(firstCode(condition.clinicalStatus) ?? ''))
    .filter((condition) => (conditionOnset(condition)?.slice(0, 10) ?? '') > day)
    .sort((left, right) => (conditionOnset(left) ?? '').localeCompare(conditionOnset(right) ?? ''));
}

/**
 * An AJCC stage or TNM category in one wording, without the word "Stage". SNOMED words them two ways
 * ("IIA (AJCC)", "American Joint Committee on Cancer stage IIA"; "cT3 (AJCC)", "American Joint Committee
 * on Cancer cT3"); both become "IIA (AJCC)" or "cT3 (AJCC)".
 */
export function stageValueLabel(value: string) {
  const ajcc = value.match(/^(?:American Joint Committee on Cancer|AJCC)\s+(?:stage\s+)?(\S+)$/i);
  return ajcc ? `${ajcc[1]} (AJCC)` : value.replace(/^stage\s+/i, '');
}

/**
 * The latest stage, with the year it was assigned so it isn't read as the current extent of disease.
 */
export function getLatestStage(resources: FhirResource[], condition: Condition | undefined): BannerItem | undefined {
  const stage = finalObservations(resources, 'mcode-tnm-stage-group', STAGE_GROUP_CODES)
    .find((observation) => observation.valueCodeableConcept && isAboutCondition(observation, condition));
  const value = conceptLabel(stage?.valueCodeableConcept);
  if (!stage || !value) {
    return undefined;
  }
  const stageDate = effectiveDate(stage);
  return {
    caption: 'Stage',
    value: stageValueLabel(value),
    note: stageDate?.slice(0, 4),
    tooltip: [conceptLabel(stage.code), formatDate(stageDate)].filter(Boolean).join(' · '),
  };
}

export function getLatestEcog(resources: FhirResource[]): BannerItem | undefined {
  const ecog = finalObservations(resources, 'mcode-ecog-performance-status', ECOG_CODES)[0];
  const value = ecog?.valueInteger ?? ecog?.valueQuantity?.value ?? conceptLabel(ecog?.valueCodeableConcept);
  if (!ecog || value === undefined) {
    return undefined;
  }
  return { caption: 'ECOG', value: String(value), tooltip: formatDate(effectiveDate(ecog)) };
}

/**
 * The disease status chip: the latest CancerDiseaseStatus about the cancer, else an informative
 * clinicalStatus of the cancer itself, else nothing. Status is never inferred from other data.
 */
export function getDiseaseStatus(resources: FhirResource[], condition: Condition | undefined): DiseaseStatus | undefined {
  const status = finalObservations(resources, 'mcode-cancer-disease-status', DISEASE_STATUS_CODES)
    .find((observation) => observation.valueCodeableConcept && isAboutCondition(observation, condition));
  if (status) {
    const code = status.valueCodeableConcept?.coding?.find((coding) => coding.system === SNOMED)?.code;
    const known = code ? DISEASE_STATUS_VALUES[code] : undefined;
    const evidence = status.extension
      ?.filter((extension) => EVIDENCE_TYPE_EXTENSIONS.includes(extension.url))
      .map((extension) => conceptLabel(extension.valueCodeableConcept))
      .filter(Boolean);
    const date = formatDate(effectiveDate(status));
    return {
      label: known?.label ?? conceptLabel(status.valueCodeableConcept) ?? 'Disease status',
      trend: known?.trend ?? 'other',
      date,
      tooltip: [
        `Cancer disease status: ${conceptLabel(status.valueCodeableConcept)}`,
        date && `Asserted ${date}`,
        evidence?.length ? `Evidence: ${evidence.join(', ')}` : undefined,
      ].filter(Boolean).join('\n'),
    };
  }

  const clinicalStatus = INFORMATIVE_CLINICAL_STATUSES[firstCode(condition?.clinicalStatus) ?? ''];
  if (clinicalStatus) {
    return { label: clinicalStatus, trend: 'other', tooltip: 'Clinical status of the primary cancer condition' };
  }
  return undefined;
}
