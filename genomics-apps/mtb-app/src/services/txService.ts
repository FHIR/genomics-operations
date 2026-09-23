// THERAPEUTIC IMPLICATION API CALL
import { FhirResponse, FhirObservation } from '@/utils/fhirInterfaces';
import { TxComponentCodes } from '@/utils/ComponentCodes';
import { processFhirResponse } from '@/utils/fhirProcessor';

// URLs
const FHIR_BASE_URL = 'https://fhir-gen-ops.herokuapp.com';
const SUBJECT_OPS_URL = `${FHIR_BASE_URL}/subject-operations/phenotype-operations/$find-subject-tx-implications`;

const CIVIC_IDENTIFIER_URLS: Record<string, (value: string) => string> = {
    'https://civicdb.org/variant': (value) => `https://civicdb.org/variants/${value}/summary`,
    'https://civicdb.org/evidence': (value) => `https://civicdb.org/evidence/${value}`,
    'https://civicdb.org/molecular-profiles': (value) => `https://civicdb.org/molecular-profiles/${value}`,
};

export interface ProcessedTxImplication {
    resourceId?: string;
    phenotypicContext: string;
    evidenceLevel: string;
    medication: string;
    therapeuticImplication: string;
    hyperlink: string;
    hyperlinkLabel?: string;
    hyperlinkTitle?: string;
    clinicalTrialId?: string;
    therapeuticImplicationDisplay?: string;
    sourceObservationIds?: string[];
}

const getTxImplicationIdentityKey = (implication: ProcessedTxImplication) => [
    implication.resourceId,
    implication.evidenceLevel,
    implication.medication,
    implication.phenotypicContext,
    implication.therapeuticImplication,
    implication.therapeuticImplicationDisplay,
    implication.clinicalTrialId,
].join('|');

const mergeSourceObservationIds = (
    primary?: string[],
    secondary?: string[]
) => {
    const mergedIds = [...(primary ?? []), ...(secondary ?? [])].filter(Boolean);

    return mergedIds.length > 0 ? [...new Set(mergedIds)].sort() : undefined;
};

const preferLongerValue = (primary?: string, secondary?: string) => {
    if (!primary) {
        return secondary;
    }

    if (!secondary) {
        return primary;
    }

    return secondary.length > primary.length ? secondary : primary;
};

const mergeTxImplication = (
    primary: ProcessedTxImplication,
    secondary: ProcessedTxImplication
): ProcessedTxImplication => ({
    ...primary,
    ...secondary,
    phenotypicContext: preferLongerValue(primary.phenotypicContext, secondary.phenotypicContext) ?? '',
    evidenceLevel: preferLongerValue(primary.evidenceLevel, secondary.evidenceLevel) ?? '',
    medication: preferLongerValue(primary.medication, secondary.medication) ?? '',
    therapeuticImplication: preferLongerValue(primary.therapeuticImplication, secondary.therapeuticImplication) ?? '',
    hyperlink: secondary.hyperlink || primary.hyperlink,
    hyperlinkLabel: secondary.hyperlinkLabel || primary.hyperlinkLabel,
    hyperlinkTitle: secondary.hyperlinkTitle || primary.hyperlinkTitle,
    clinicalTrialId: secondary.clinicalTrialId || primary.clinicalTrialId,
    therapeuticImplicationDisplay: secondary.therapeuticImplicationDisplay || primary.therapeuticImplicationDisplay,
    sourceObservationIds: mergeSourceObservationIds(primary.sourceObservationIds, secondary.sourceObservationIds),
});

export const getTxImplicationDedupKey = (implication: ProcessedTxImplication) => [
    getTxImplicationIdentityKey(implication),
    (implication.sourceObservationIds ?? []).slice().sort().join(','),
].join('|');

export const dedupeTxImplications = (implications: ProcessedTxImplication[]) => {
    const byIdentity = new Map<string, ProcessedTxImplication>();

    implications.forEach((implication) => {
        const key = getTxImplicationIdentityKey(implication);
        const existing = byIdentity.get(key);

        if (!existing) {
            byIdentity.set(key, implication);
            return;
        }

        byIdentity.set(key, mergeTxImplication(existing, implication));
    });

    return Array.from(byIdentity.values());
};

const findComponent = (resource: FhirObservation, code: string) =>
    resource.component?.find(component =>
        component.code?.coding?.some(coding => coding.code === code)
    );

const getPrimaryCoding = (resource: FhirObservation, code: string) =>
    findComponent(resource, code)?.valueCodeableConcept?.coding?.[0];

const getEvidenceLevelValue = (resource: FhirObservation) => {
    const evidenceComponent = findComponent(resource, TxComponentCodes.evidenceLevel);

    return evidenceComponent?.valueCodeableConcept?.text
        || evidenceComponent?.valueCodeableConcept?.coding?.[0]?.display
        || evidenceComponent?.valueCodeableConcept?.coding?.[0]?.code
        || '';
};

const getDerivedFromObservationIds = (resource: FhirObservation) =>
    (resource.derivedFrom ?? [])
        .map(reference => reference.reference?.match(/^Observation\/(.+)$/)?.[1] ?? '')
        .filter(Boolean);

const getCivicHyperlink = (resource: FhirObservation) => {
    const preferredSystems = [
        'https://civicdb.org/evidence',
        'https://civicdb.org/variant',
        'https://civicdb.org/molecular-profiles',
    ];

    const identifier = preferredSystems
        .map((system) => resource.identifier?.find((candidate) => candidate.system === system))
        .find(Boolean);

    if (!identifier?.system || !identifier.value) {
        return '';
    }

    const buildUrl = CIVIC_IDENTIFIER_URLS[identifier.system];
    return buildUrl ? buildUrl(identifier.value) : '';
};

// Helper to extract CIViC hyperlink
const txExtraFieldsExtractor = (resource: FhirObservation): Partial<ProcessedTxImplication> => {
    const evidenceLevel = getEvidenceLevelValue(resource);
    const therapeuticImplicationCoding = getPrimaryCoding(resource, TxComponentCodes.therapeuticImplication);
    const medicationCoding = getPrimaryCoding(resource, TxComponentCodes.medication);

    if (/^clinical trial$/i.test(evidenceLevel)) {
        const clinicalTrialId = medicationCoding?.code || '';

        return {
            resourceId: resource.id,
            hyperlink: clinicalTrialId ? `https://clinicaltrials.gov/study/${clinicalTrialId}` : '',
            hyperlinkLabel: clinicalTrialId ? 'Clinical Trial link' : undefined,
            hyperlinkTitle: medicationCoding?.display || '',
            clinicalTrialId,
            therapeuticImplicationDisplay: therapeuticImplicationCoding?.display || therapeuticImplicationCoding?.code || '',
            sourceObservationIds: getDerivedFromObservationIds(resource),
        };
    }

    return {
        resourceId: resource.id,
        hyperlink: getCivicHyperlink(resource),
        sourceObservationIds: getDerivedFromObservationIds(resource),
    };
};

export async function processTxImplications(
    variant: string,
    subjectId: string,
    experimental = false
): Promise<ProcessedTxImplication[]> {
    try {
        const response = await fetch(
            `${SUBJECT_OPS_URL}?subject=${encodeURIComponent(subjectId)}&variants=${encodeURIComponent(variant)}&experimental=${experimental}`
        );

        if (!response.ok) {
            return [];
        }

        const data: FhirResponse = await response.json();
        return dedupeTxImplications(
            processFhirResponse<ProcessedTxImplication>(
                data,
                TxComponentCodes,
                'https://civicdb.org/variant',
                txExtraFieldsExtractor
            )
        );
    } catch (error) {
        console.error(`Error processTxImplications:`, error);
        return [];
    }
}

export async function processTxImplicationsForRange(
    range: string,
    subjectId: string,
    experimental = false
): Promise<ProcessedTxImplication[]> {
    try {
        const response = await fetch(
            `${SUBJECT_OPS_URL}?subject=${encodeURIComponent(subjectId)}&ranges=${encodeURIComponent(range)}&experimental=${experimental}`
        );

        if (!response.ok) {
            return [];
        }

        const data: FhirResponse = await response.json();
        return dedupeTxImplications(
            processFhirResponse<ProcessedTxImplication>(
                data,
                TxComponentCodes,
                'https://civicdb.org/variant',
                txExtraFieldsExtractor
            )
        );
    } catch (error) {
        console.error(`Error processTxImplicationsForRange:`, error);
        return [];
    }
}
