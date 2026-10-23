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
  [key: string]: unknown;
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
  onsetDateTime?: string;
  recordedDate?: string;
}

export interface Observation extends FhirResource {
  resourceType: 'Observation';
  status?: string;
  category?: CodeableConcept[];
  code?: CodeableConcept;
  focus?: Reference[];
  effectiveDateTime?: string;
  effectivePeriod?: Period;
  issued?: string;
  valueCodeableConcept?: CodeableConcept;
  valueInteger?: number;
  valueQuantity?: Quantity;
  valueString?: string;
}

export interface AllergyIntolerance extends FhirResource {
  resourceType: 'AllergyIntolerance';
  clinicalStatus?: CodeableConcept;
  verificationStatus?: CodeableConcept;
  code?: CodeableConcept;
}

export interface Bundle {
  resourceType: 'Bundle';
  total?: number;
  link?: { relation: string; url: string }[];
  entry?: { resource?: FhirResource }[];
}
