import { Bundle, FhirResource, Patient } from '@/types/fhir';

// Clinical data (mCODE and anything else) lives on a separate FHIR server from the genomic data
export const CLINICAL_FHIR_BASE_URL = 'https://hapi.fhir.org/baseR4';
export const MRN_IDENTIFIER_SYSTEM = 'http://elimu.io/fhir/mrn';

export interface ClinicalRecord {
  mrn: string;
  patient: Patient;
  resources: FhirResource[];
}

// Cached per MRN for the session; a failed request is evicted so it can be retried
const recordCache = new Map<string, Promise<ClinicalRecord | null>>();
let clinicalMrnsPromise: Promise<string[]> | null = null;

async function fetchBundle(url: string): Promise<Bundle> {
  const response = await fetch(url, { headers: { Accept: 'application/fhir+json' } });
  if (!response.ok) {
    throw new Error(`FHIR request failed (HTTP ${response.status}): ${url}`);
  }
  return response.json();
}

// Fetch a search or operation result, following `next` links until all pages are read
async function fetchAllPages(url: string): Promise<FhirResource[]> {
  const resources: FhirResource[] = [];
  let nextUrl: string | undefined = url;
  while (nextUrl) {
    const bundle = await fetchBundle(nextUrl);
    for (const entry of bundle.entry ?? []) {
      if (entry.resource) {
        resources.push(entry.resource);
      }
    }
    nextUrl = bundle.link?.find((link) => link.relation === 'next')?.url;
  }
  return resources;
}

async function loadClinicalRecord(mrn: string): Promise<ClinicalRecord | null> {
  const identifier = encodeURIComponent(`${MRN_IDENTIFIER_SYSTEM}|${mrn}`);
  const patients = (await fetchAllPages(`${CLINICAL_FHIR_BASE_URL}/Patient?identifier=${identifier}`))
    .filter((resource): resource is Patient => resource.resourceType === 'Patient');
  if (patients.length === 0) {
    return null;
  }
  if (patients.length > 1) {
    console.warn(`${patients.length} patients on the FHIR server have MRN ${mrn}; using Patient/${patients[0].id}`);
  }

  const patient = patients[0];
  const everything = await fetchAllPages(`${CLINICAL_FHIR_BASE_URL}/Patient/${patient.id}/$everything`);
  // $everything can return the same resource on more than one page
  const unique = new Map(everything.map((resource) => [`${resource.resourceType}/${resource.id}`, resource]));
  return { mrn, patient, resources: [...unique.values()] };
}

/**
 * The patient's whole clinical record, or null when the FHIR server has no patient with this MRN.
 */
export function fetchClinicalRecord(mrn: string): Promise<ClinicalRecord | null> {
  let promise = recordCache.get(mrn);
  if (!promise) {
    promise = loadClinicalRecord(mrn);
    recordCache.set(mrn, promise);
    promise.catch(() => recordCache.delete(mrn));
  }
  return promise;
}

/**
 * Every MRN that has a patient on the clinical FHIR server, detected once per page load.
 */
export function fetchClinicalMrns(): Promise<string[]> {
  if (!clinicalMrnsPromise) {
    const identifier = encodeURIComponent(`${MRN_IDENTIFIER_SYSTEM}|`);
    clinicalMrnsPromise = fetchAllPages(`${CLINICAL_FHIR_BASE_URL}/Patient?identifier=${identifier}&_elements=identifier&_count=1000`)
      .then((resources) => Array.from(new Set(
        (resources as Patient[]).flatMap((patient) => (patient.identifier ?? [])
          .filter((id) => id.system === MRN_IDENTIFIER_SYSTEM && id.value)
          .map((id) => id.value as string))
      )));
    clinicalMrnsPromise.catch(() => {
      clinicalMrnsPromise = null;
    });
  }
  return clinicalMrnsPromise;
}
