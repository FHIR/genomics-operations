// Minimal FHIR R4 shapes for the clinical data the app reads. Only the elements the app uses are typed;
// everything else stays reachable through the index signature.

export interface Coding {
  system?: string;
  code?: string;
  display?: string;
}

export interface CodeableConcept {
  coding?: Coding[];
  text?: string;
}

export interface Reference {
  reference?: string;
  display?: string;
}

export interface Period {
  start?: string;
  end?: string;
}

export interface Quantity {
  value?: number;
  unit?: string;
  code?: string;
}

export interface Extension {
  url: string;
  valueCodeableConcept?: CodeableConcept;
  valueReference?: Reference;
  [key: string]: unknown;
}

export interface Range {
  low?: Quantity;
  high?: Quantity;
}

// A coded element that may carry extensions, e.g. Condition.bodySite with a laterality qualifier
export interface ExtendedCodeableConcept extends CodeableConcept {
  extension?: Extension[];
}

// The value[x] choices the app displays
export interface ObservationValue {
  valueCodeableConcept?: CodeableConcept;
  valueInteger?: number;
  valueQuantity?: Quantity;
  valueString?: string;
  valueBoolean?: boolean;
  valueRange?: Range;
  valueRatio?: { numerator?: Quantity; denominator?: Quantity };
  valueDateTime?: string;
  dataAbsentReason?: CodeableConcept;
}

export interface FhirResource {
  resourceType: string;
  id?: string;
  meta?: { profile?: string[]; lastUpdated?: string };
  extension?: Extension[];
  [key: string]: unknown;
}

export interface Identifier {
  system?: string;
  value?: string;
}

export interface Patient extends FhirResource {
  resourceType: 'Patient';
  identifier?: Identifier[];
  name?: { use?: string; text?: string; family?: string; given?: string[] }[];
  gender?: string;
  birthDate?: string;
}

export interface Condition extends FhirResource {
  resourceType: 'Condition';
  clinicalStatus?: CodeableConcept;
  verificationStatus?: CodeableConcept;
  category?: CodeableConcept[];
  code?: CodeableConcept;
  bodySite?: ExtendedCodeableConcept[];
  onsetDateTime?: string;
  onsetPeriod?: Period;
  abatementDateTime?: string;
  recordedDate?: string;
  stage?: { summary?: CodeableConcept; type?: CodeableConcept; assessment?: Reference[] }[];
}

export interface Observation extends FhirResource, ObservationValue {
  resourceType: 'Observation';
  status?: string;
  category?: CodeableConcept[];
  code?: CodeableConcept;
  focus?: Reference[];
  effectiveDateTime?: string;
  effectivePeriod?: Period;
  issued?: string;
  method?: CodeableConcept;
  interpretation?: CodeableConcept[];
  referenceRange?: { low?: Quantity; high?: Quantity; text?: string }[];
  hasMember?: Reference[];
  component?: (ObservationValue & { code?: CodeableConcept })[];
}

export interface AllergyIntolerance extends FhirResource {
  resourceType: 'AllergyIntolerance';
  clinicalStatus?: CodeableConcept;
  verificationStatus?: CodeableConcept;
  category?: string[];
  criticality?: string;
  code?: CodeableConcept;
  reaction?: { manifestation?: CodeableConcept[]; severity?: string }[];
}

export interface DiagnosticReport extends FhirResource {
  resourceType: 'DiagnosticReport';
  status?: string;
  category?: CodeableConcept[];
  code?: CodeableConcept;
  effectiveDateTime?: string;
  effectivePeriod?: Period;
  issued?: string;
  result?: Reference[];
  conclusion?: string;
  conclusionCode?: CodeableConcept[];
  presentedForm?: { contentType?: string; data?: string; title?: string }[];
}

export interface FamilyMemberHistory extends FhirResource {
  resourceType: 'FamilyMemberHistory';
  status?: string;
  date?: string;
  relationship?: CodeableConcept;
  deceasedBoolean?: boolean;
  condition?: {
    code?: CodeableConcept;
    onsetAge?: Quantity;
    onsetString?: string;
    contributedToDeath?: boolean;
  }[];
}

export interface Bundle {
  resourceType: 'Bundle';
  total?: number;
  link?: { relation: string; url: string }[];
  entry?: { resource?: FhirResource }[];
}
