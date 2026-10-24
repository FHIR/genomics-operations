// The case synthesis for a patient. For the demo it is pre-generated Markdown stored in
// public/data/synthesis/<MRN>.md (mtb-clinical-data-design.md section 6.2); a real AI call (with a
// server-held key) could replace this function without changing the Synthesis tab.

const SYNTHESIS_BASE_PATH = '/data/synthesis';

/**
 * The synthesis as Markdown, or null when there is none for this patient.
 */
export async function fetchSynthesis(mrn: string): Promise<string | null> {
  const response = await fetch(`${SYNTHESIS_BASE_PATH}/${encodeURIComponent(mrn)}.md`);
  if (response.status === 404) {
    return null;
  }
  if (!response.ok) {
    throw new Error(`Synthesis request failed (HTTP ${response.status})`);
  }
  return response.text();
}
