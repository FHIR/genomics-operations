export type FhirCoding = {
    display?: string;
    text?: string;
    code?: string;
};

export type FhirCodeableConcept = {
    coding?: FhirCoding[];
    text?: string;
};

export type FhirAnnotation = {
    text?: string;
};

export type FhirObservationComponent = {
    code?: {
        text?: string;
    };
    valueInteger?: number;
    interpretation?: FhirCodeableConcept | FhirCodeableConcept[];
};

export type FhirObservation = {
    resourceType?: string;
    valueInteger?: number;
    interpretation?: FhirCodeableConcept | FhirCodeableConcept[];
    note?: FhirAnnotation[];
    component?: FhirObservationComponent[];
    issued?: string;
    dataAbsentReason?: FhirCodeableConcept;
};

export type OncogenicityGaugeKind =
    | 'dark-green'
    | 'light-green'
    | 'yellow'
    | 'orange'
    | 'red'
    | 'undetermined';

export type OncogenicityEvidenceStatus = 'idle' | 'loading' | 'ready' | 'error';

export interface OncogenicityPredictionResult {
    key: string;
    spdi: string;
    hgvs?: string;
    status: 'loading' | 'ready';
    gauge: OncogenicityGaugeKind;
    score?: number;
    interpretation?: string;
    observation?: FhirObservation;
    errorMessage?: string;
    evidenceStatus: OncogenicityEvidenceStatus;
    evidenceJson?: unknown;
    evidenceError?: string;
}

function getFirstConcept(concept?: FhirCodeableConcept | FhirCodeableConcept[]): FhirCodeableConcept | undefined {
    if (!concept) {
        return undefined;
    }

    return Array.isArray(concept) ? concept[0] : concept;
}

export function getConceptDisplayText(concept?: FhirCodeableConcept | FhirCodeableConcept[]) {
    const resolvedConcept = getFirstConcept(concept);
    if (!resolvedConcept) {
        return undefined;
    }

    const firstCoding = resolvedConcept.coding?.[0];
    return firstCoding?.display || firstCoding?.text || resolvedConcept.text;
}

export function getConceptNarrativeText(concept?: FhirCodeableConcept | FhirCodeableConcept[]) {
    const resolvedConcept = getFirstConcept(concept);
    if (!resolvedConcept) {
        return undefined;
    }

    const firstCoding = resolvedConcept.coding?.[0];
    return resolvedConcept.text || firstCoding?.text || firstCoding?.display;
}

export function getObservationCaveat(observation?: FhirObservation) {
    return observation?.note?.map((note) => note.text).filter(Boolean).join('\n');
}

export function getGaugeKindForScore(score?: number): OncogenicityGaugeKind {
    if (typeof score !== 'number' || Number.isNaN(score)) {
        return 'undetermined';
    }

    if (score <= -7) {
        return 'dark-green';
    }

    if (score >= -6 && score <= -1) {
        return 'light-green';
    }

    if (score >= 0 && score <= 5) {
        return 'yellow';
    }

    if (score >= 6 && score <= 9) {
        return 'orange';
    }

    return 'red';
}

export function getGaugeImagePath(gauge: OncogenicityGaugeKind) {
    switch (gauge) {
        case 'dark-green':
            return '/gauge_dark_green.png';
        case 'light-green':
            return '/gauge_light_green.png';
        case 'yellow':
            return '/gauge_yellow.png';
        case 'orange':
            return '/gauge_orange.png';
        case 'red':
            return '/gauge_red.png';
        case 'undetermined':
        default:
            return '/gauge_undetermined.png';
    }
}
