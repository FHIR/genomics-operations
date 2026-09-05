import { processFhirResponse } from '@/utils/fhirProcessor';
import { dxComponentCodes } from '@/utils/ComponentCodes';
import { FhirResponse, FhirObservation } from '@/utils/fhirInterfaces';

// URLs
const FHIR_BASE_URL = 'https://fhir-gen-ops.herokuapp.com';
const SUBJECT_OPS_URL = `${FHIR_BASE_URL}/subject-operations/phenotype-operations/$find-subject-dx-implications`;

const CLINVAR_SYSTEM = 'http://www.ncbi.nlm.nih.gov/clinvar';

export interface DxImplication {
    clinicalSignificance?: string;
    predictedPhenotype?: string;
    evidenceLevel?: string;
    variantId?: string;
    clinvarLink?: string;
    stars?: number | null;
    starsEmoji?: string | null;
}

// Helper to extract non-component fields (e.g., variantId & link)
const dxExtraFieldsExtractor = (resource: FhirObservation): Partial<DxImplication> => {
    const variantId = resource.identifier?.find(id => id.system === CLINVAR_SYSTEM)?.value || '';
    const clinvarLink = variantId ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${variantId}/` : '';
    return { variantId, clinvarLink };
};

export async function processDxImplications(
    variant: string,
    subjectId: string
): Promise<DxImplication[]> {
    try {
        const response = await fetch(
            `${SUBJECT_OPS_URL}?subject=${encodeURIComponent(subjectId)}&variants=${encodeURIComponent(variant)}`
        );

        if (!response.ok) {
            return [];
        }

        const data: FhirResponse = await response.json();

        return await processFhirResponse<DxImplication>(
            data,
            dxComponentCodes,
            CLINVAR_SYSTEM,
            dxExtraFieldsExtractor
        );
    } catch (error) {
        console.error(`Error processDxImplications:`, error);
        return [];
    }
}