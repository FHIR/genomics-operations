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

export async function fetchPredictionObservation(variant: string, tumorType?: string) {
    const params = new URLSearchParams({ variant });
    if (tumorType) {
        params.set('tumorType', tumorType);
    }
    const response = await fetchJsonResponse(`${PREDICTOR_BASE_URL}/predictOncogenicity?${params.toString()}`);

    if (!response.ok) {
        throw new Error(`Request failed: ${response.status} ${response.statusText}`);
    }

    return response.json() as Promise<unknown>;
}

export async function fetchEvidenceSummary(variant: string, tumorType?: string) {
    const params = new URLSearchParams({ variant });
    if (tumorType) {
        params.set('tumorType', tumorType);
    }
    const response = await fetchJsonResponse(`${PREDICTOR_BASE_URL}/summarizeEvidence?${params.toString()}`);

    if (!response.ok) {
        throw new Error(`Request failed: ${response.status} ${response.statusText}`);
    }

    return response.json() as Promise<unknown>;
}
