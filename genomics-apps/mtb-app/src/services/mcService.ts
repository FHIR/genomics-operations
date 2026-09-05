// MOLECULAR CONSEQUENCE API CALL
import { FhirResponse, FhirObservation } from '@/utils/fhirInterfaces';
import { McComponentCodes } from '@/utils/ComponentCodes';
import { processFhirResponse } from '@/utils/fhirProcessor';

// URLs
const FHIR_BASE_URL = 'https://fhir-gen-ops.herokuapp.com';
const SUBJECT_OPS_URL = `${FHIR_BASE_URL}/subject-operations/phenotype-operations/$find-subject-molecular-consequences`;

export interface MolecularConsequence {
    featureConsequence: string;
    interpretation: string;
    impact: string;
    hyperlink: string;
    proteinChange?: string;
}

// Helper to extract interpretation and add fixed hyperlink
const mcExtraFieldsExtractor = (resource: FhirObservation): Partial<MolecularConsequence> => {
    // Get the impact from the interpretation text
    const impact = resource.interpretation?.[0]?.text?.toUpperCase() || 'UNKNOWN';

    return {
        interpretation: impact,
        impact,
        hyperlink: 'https://pcingola.github.io/SnpEff/'
    };
};

export async function processMolecularConsequences(
    variant: string,
    subjectId: string
): Promise<MolecularConsequence[]> {
    try {
        const response = await fetch(
            `${SUBJECT_OPS_URL}?subject=${encodeURIComponent(subjectId)}&variants=${encodeURIComponent(variant)}`
        );

        if (!response.ok) {
            return [];
        }

        const data: FhirResponse = await response.json();
        return processFhirResponse<MolecularConsequence>(
            data,
            McComponentCodes,
            undefined,
            mcExtraFieldsExtractor
        );
    } catch (error) {
        console.error(`Error processMolecularConsequences:`, error);
        return [];
    }
}
