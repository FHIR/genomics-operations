// Slide studies usually hold several slides, and IDC's slide viewer opens one of its own choosing, often a frozen
// section. The FFPE diagnostic slide is the one pathologists read, and the one IDC's AI overlays (nuclei, TIL maps)
// are drawn on. The record lists each slide's series UID but not its name, so the name comes from IDC.

const SERIES_DESCRIPTION_TAG = '0008103E';
const SERIES_UID_TAG = '0020000E';

type DicomJson = Record<string, { Value?: unknown[] }>;

// Cached per study for the session; a failed lookup is evicted so a later click can retry
const cache = new Map<string, Promise<string | undefined>>();

async function lookUp(dicomWebBase: string, studyUid: string): Promise<string | undefined> {
  const url = `${dicomWebBase}/studies/${encodeURIComponent(studyUid)}/series?Modality=SM&includefield=${SERIES_DESCRIPTION_TAG}`;
  const response = await fetch(url, { headers: { Accept: 'application/dicom+json' } });
  if (!response.ok) {
    throw new Error(`IDC series lookup failed (HTTP ${response.status})`);
  }
  const series = ((await response.json()) as DicomJson[]).map((item) => ({
    uid: item[SERIES_UID_TAG]?.Value?.[0] as string | undefined,
    description: String(item[SERIES_DESCRIPTION_TAG]?.Value?.[0] ?? ''),
  }));
  // TCGA names its slides e.g. "FFPE HE TP DX1" (diagnostic) or "Frozen HE TP TS1"
  const ffpe = series.filter((item) => /\bFFPE\b/i.test(item.description));
  return (ffpe.find((item) => /\bDX\d*\b/i.test(item.description)) ?? ffpe[0])?.uid;
}

/**
 * The series UID of the study's FFPE diagnostic slide, or undefined when there isn't one.
 */
export function findDiagnosticSlide(dicomWebBase: string, studyUid: string): Promise<string | undefined> {
  const key = `${dicomWebBase}|${studyUid}`;
  let promise = cache.get(key);
  if (!promise) {
    promise = lookUp(dicomWebBase, studyUid);
    cache.set(key, promise);
    promise.catch(() => cache.delete(key));
  }
  return promise;
}
