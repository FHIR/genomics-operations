// THERAPEUTIC IMPLICATION API CALL
import { FhirResponse, FhirObservation } from '@/utils/fhirInterfaces';
import { TxComponentCodes } from '@/utils/ComponentCodes';
import { processFhirResponse } from '@/utils/fhirProcessor';

// URLs
const FHIR_BASE_URL = 'https://fhir-gen-ops.herokuapp.com';
const SUBJECT_OPS_URL = `${FHIR_BASE_URL}/subject-operations/phenotype-operations/$find-subject-tx-implications`;

export interface ProcessedTxImplication {
    phenotypicContext: string;
    evidenceLevel: string;
    medication: string;
    therapeuticImplication: string;
    hyperlink: string;
}

// Helper to extract CIViC hyperlink
const txExtraFieldsExtractor = (resource: FhirObservation): Partial<ProcessedTxImplication> => {
    const civicIdentifier = resource.identifier?.find(
        id => id.system === 'https://civicdb.org/variant'
    );

    const variantId = civicIdentifier?.value;
    const url = variantId ? `https://civicdb.org/variants/${variantId}/summary` : '';

    return {
        hyperlink: url
    };
};

export async function processTxImplications(
    variant: string,
    subjectId: string
): Promise<ProcessedTxImplication[]> {
    try {
        const response = await fetch(
            `${SUBJECT_OPS_URL}?subject=${encodeURIComponent(subjectId)}&variants=${encodeURIComponent(variant)}`
        );

        if (!response.ok) {
            return [];
        }

        const data: FhirResponse = await response.json();
        return processFhirResponse<ProcessedTxImplication>(
            data,
            TxComponentCodes,
            'https://civicdb.org/variant',
            txExtraFieldsExtractor
        );
    } catch (error) {
        console.error(`Error processTxImplications:`, error);
        return [];
    }
}
