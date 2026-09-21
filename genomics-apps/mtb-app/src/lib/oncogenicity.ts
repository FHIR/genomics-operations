const NCBI_VARIATION_BASE_URL = 'https://api.ncbi.nlm.nih.gov/variation/v0';
const PREDICTOR_BASE_URL = 'https://oncogenicity-predictor.onrender.com';

async function fetchJsonResponse(url: string) {
    const response = await fetch(url, {
        method: 'GET',
        headers: {
            Accept: 'application/json',
        },
        cache: 'no-store',
    });

    return response;
}

async function fetchJson<T>(url: string): Promise<T> {
    const response = await fetchJsonResponse(url);

    if (!response.ok) {
        throw new Error(`Request failed: ${response.status} ${response.statusText}`);
    }

    return response.json() as Promise<T>;
}

export async function convertSpdiToHgvs(spdi: string) {
    const url = `${NCBI_VARIATION_BASE_URL}/spdi/${encodeURIComponent(spdi)}/hgvs`;
    const response = await fetchJson<{ data?: { hgvs?: string } }>(url);
    const hgvs = response.data?.hgvs;

    if (!hgvs) {
        throw new Error('HGVS conversion did not return a value');
    }

    return hgvs;
}

export async function fetchPredictionObservation(hgvs: string, tumorType?: string) {
    const params = new URLSearchParams({ variant: hgvs });
    if (tumorType) {
        params.set('tumorType', tumorType);
    }
    const response = await fetchJsonResponse(`${PREDICTOR_BASE_URL}/predictOncogenicity?${params.toString()}`);

    if (!response.ok) {
        throw new Error(`Request failed: ${response.status} ${response.statusText}`);
    }

    return response.json() as Promise<unknown>;
}

export async function fetchEvidenceSummary(hgvs: string, tumorType?: string) {
    const params = new URLSearchParams({ variant: hgvs });
    if (tumorType) {
        params.set('tumorType', tumorType);
    }
    const response = await fetchJsonResponse(`${PREDICTOR_BASE_URL}/summarizeEvidence?${params.toString()}`);

    if (!response.ok) {
        throw new Error(`Request failed: ${response.status} ${response.statusText}`);
    }

    return response.json() as Promise<unknown>;
}
