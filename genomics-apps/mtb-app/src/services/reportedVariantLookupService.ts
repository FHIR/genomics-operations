import { getFeatureCoordinates } from './variantService';

// Lookups for layering reported variants onto the Genomics Explorer: a reported variant's SPDI in each build,
// and a gene's or transcript's range. They're made only for the few reported variants that need them, once each, and cached for the session.

const UTILITIES_URL = 'https://fhir-gen-ops.herokuapp.com/utilities';
const GENOMIC_RANGE = /^NC_\d+\.\d+:\d+-\d+$/;

export interface NormalizedVariant {
  b37SPDI?: string;
  b38SPDI?: string;
}

const normalizedCache = new Map<string, Promise<NormalizedVariant | null>>();
const geneCache = new Map<string, Promise<string | undefined>>();
const rangeCache = new Map<string, Promise<string | undefined>>();

// A failed lookup is evicted so a later search can retry it
function cached<T>(cache: Map<string, Promise<T>>, key: string, load: () => Promise<T>, fallback: T): Promise<T> {
  let promise = cache.get(key);
  if (!promise) {
    promise = load().catch((error) => {
      console.error(`Lookup failed for ${key}:`, error);
      cache.delete(key);
      return fallback;
    });
    cache.set(key, promise);
  }
  return promise;
}

/**
 * Build 37 and build 38 SPDIs for a variant given as SPDI, gHGVS or c.HGVS (the genomics server's
 * normalize-variant, which uses NCBI Variation Services).
 */
export function normalizeVariant(variant: string): Promise<NormalizedVariant | null> {
  return cached(normalizedCache, variant, async () => {
    const response = await fetch(`${UTILITIES_URL}/normalize-variant?variant=${encodeURIComponent(variant)}`);
    if (!response.ok) {
      throw new Error(`normalize-variant failed (HTTP ${response.status})`);
    }
    const data = await response.json() as NormalizedVariant;
    return { b37SPDI: data.b37SPDI, b38SPDI: data.b38SPDI };
  }, null);
}

/**
 * The gene symbol of a transcript, e.g. NM_000455.5 is STK11 (version numbers are ignored).
 */
export function geneForTranscript(transcript: string): Promise<string | undefined> {
  return cached(geneCache, transcript, async () => {
    const response = await fetch(`${UTILITIES_URL}/get-feature-coordinates?transcript=${encodeURIComponent(transcript)}`);
    if (!response.ok) {
      throw new Error(`get-feature-coordinates failed (HTTP ${response.status})`);
    }
    const data = await response.json() as { geneSymbol?: string }[];
    return data?.[0]?.geneSymbol;
  }, undefined);
}

/**
 * The genomic range a search term covers: a range as typed, or a gene's coordinates (as the search itself resolves them).
 */
export function resolveSearchRange(term: string): Promise<string | undefined> {
  if (GENOMIC_RANGE.test(term)) {
    return Promise.resolve(term);
  }
  return cached(rangeCache, term, () => getFeatureCoordinates(term), undefined);
}
