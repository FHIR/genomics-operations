import { StructuralDetails, Variant } from '@/types/variants';
import { getStructuralDescription } from '@/lib/variantDisplay';

/**
 * The genomics server has no data for this patient (it answers 404 "Patient (…) not found").
 * The patient may still have a clinical record, and reported variants, on the clinical FHIR server.
 */
export class PatientNotOnServerError extends Error {
    constructor(subjectId: string) {
        super(`The genomics server has no data for patient ${subjectId}`);
        this.name = 'PatientNotOnServerError';
    }
}

// API Configuration
const API_CONFIG = {
    baseUrl: 'https://fhir-gen-ops.herokuapp.com',
    defaultSubjectId: 'L2345', // Example patient identifier
    endpoints: {
        utilities: '/utilities',
        subjectOps: '/subject-operations/genotype-operations'
    }
} as const;

// FHIR Resource Interfaces
interface GeneInfoResponse {
    geneId: string;
    geneSymbol: string;
    geneLink: string;
    build37Coordinates: string;
    build38Coordinates: string;
    transcripts: string[];
    MANE: string[];
}

interface FhirVariantResource {
    id?: string;
    component?: {
        code?: {
            coding?: {
                code: string;
                display?: string;
            }[];
        };
        valueCodeableConcept?: {
            coding?: {
                code: string;
                display?: string;
            }[];
            text?: string;
        };
        valueQuantity?: {
            value?: number;
        };
        valueString?: string;
        valueRange?: {
            low?: {
                value?: number;
            };
            high?: {
                value?: number;
            };
        };
    }[];
}

interface FhirParameter {
    name: string;
    part?: {
        name: string;
        resource?: FhirVariantResource;
    }[];
}

interface FhirResponse {
    parameter?: FhirParameter[];
}

// Utility Functions
const isGenomicCoordinate = (input: string): boolean => /^NC_\d+\.\d+:\d+-\d+$/.test(input);
const isGeneSymbol = (input: string): boolean => /^[A-Z0-9]+$/.test(input);

/**
 * Convert gene symbol to genomic coordinates
 */
export const getFeatureCoordinates = async (geneSymbol: string): Promise<string> => {
    try {
        // Return if already in correct format
        if (isGenomicCoordinate(geneSymbol)) {
            return geneSymbol;
        }

        const response = await fetch(
            `${API_CONFIG.baseUrl}${API_CONFIG.endpoints.utilities}/get-feature-coordinates?gene=${encodeURIComponent(geneSymbol)}`
        );

        if (!response.ok) {
            throw new Error(`API error: ${response.status} ${response.statusText}`);
        }

        const data = await response.json() as GeneInfoResponse[];

        if (!data?.[0]?.build37Coordinates) {
            throw new Error(`No genomic coordinates found for gene: ${geneSymbol}`);
        }

        return data[0].build37Coordinates;

    } catch (error) {
        console.error(`Error fetching coordinates for gene ${geneSymbol}:`, error);
        throw error;
    }
};

/**
 * Extract variant string from FHIR resource
 */
const extractVariantString = (resource: FhirVariantResource): string => {
    const variantComponent = resource.component?.find(
        c => c.code?.coding?.[0]?.code === "81252-9"
    );

    return variantComponent?.valueCodeableConcept?.coding?.[0]?.code || "Unknown variant";
};

const findComponent = (resource: FhirVariantResource, code: string) =>
    resource.component?.find(component =>
        component.code?.coding?.some(coding => coding.code === code)
    );

// A structural variant's fields: change type, location (outer range, else inner), copy number and inserted bases
const extractStructuralDetails = (resource: FhirVariantResource): StructuralDetails => {
    const outerRange = findComponent(resource, '81301-4')?.valueRange;
    const innerRange = findComponent(resource, '81302-2')?.valueRange;
    const range = outerRange?.low?.value !== undefined && outerRange?.high?.value !== undefined ? outerRange : innerRange;
    const insertedSequence = findComponent(resource, '69551-0')?.valueString;
    return {
        changeType: findComponent(resource, '48019-4')?.valueCodeableConcept?.coding?.[0]?.display,
        referenceSequence: findComponent(resource, '48013-7')?.valueCodeableConcept?.coding?.[0]?.code,
        start: range?.low?.value,
        end: range?.high?.value,
        copyNumber: findComponent(resource, '82155-3')?.valueQuantity?.value,
        // Only literal bases; a symbolic allele such as <INS> is not a sequence
        insertedSequence: insertedSequence && /^[ACGTN]+$/i.test(insertedSequence) ? insertedSequence : undefined,
    };
};

const extractGenomicSourceClass = (resource: FhirVariantResource): string =>
    findComponent(resource, '48002-0')?.valueCodeableConcept?.coding?.[0]?.display || '<unknown>';

const extractVariantAlleleFrequency = (resource: FhirVariantResource): number | undefined =>
    findComponent(resource, '81258-6')?.valueQuantity?.value;

const extractVariantResources = (data: FhirResponse): FhirVariantResource[] =>
    data.parameter
        ?.find(p => p.name === 'variants')
        ?.part
        ?.filter((p): p is { name: string; resource: FhirVariantResource } => p.name === 'variant' && !!p.resource)
        ?.map(p => p.resource) || [];

/**
 * API call to find variants for a given range
 */
export const findSubjectVariants = async (
    range: string,
    subjectId: string = API_CONFIG.defaultSubjectId
): Promise<Variant[]> => {
    try {
        if (!isGenomicCoordinate(range) && !isGeneSymbol(range)) {
            throw new Error(`Invalid input format: ${range}`);
        }

        // Convert gene symbol to genomic coordinates if needed
        const genomicRange = isGenomicCoordinate(range) ? range : await getFeatureCoordinates(range);

        const simpleVariantsUrl = `${API_CONFIG.baseUrl}${API_CONFIG.endpoints.subjectOps}/$find-subject-variants?subject=${subjectId}&ranges=${encodeURIComponent(genomicRange)}&includeVariants=true`;
        const structuralVariantsUrl = `${API_CONFIG.baseUrl}${API_CONFIG.endpoints.subjectOps}/$find-subject-structural-intersecting-variants?subject=${subjectId}&ranges=${encodeURIComponent(genomicRange)}&includeVariants=true`;

        const [simpleResponse, structuralResponse] = await Promise.all([
            fetch(simpleVariantsUrl),
            fetch(structuralVariantsUrl)
        ]);

        if (simpleResponse.status === 404 || structuralResponse.status === 404) {
            throw new PatientNotOnServerError(subjectId);
        }

        if (!simpleResponse.ok) {
            throw new Error(`API error: ${simpleResponse.status} ${simpleResponse.statusText}`);
        }

        if (!structuralResponse.ok) {
            throw new Error(`API error: ${structuralResponse.status} ${structuralResponse.statusText}`);
        }

        const [simpleData, structuralData] = await Promise.all([
            simpleResponse.json() as Promise<FhirResponse>,
            structuralResponse.json() as Promise<FhirResponse>
        ]);

        return [
            ...extractVariantResources(simpleData).map(resource => ({
                range,
                resolvedRange: genomicRange,
                sourceObservationId: resource.id,
                variant: extractVariantString(resource),
                variantType: 'simple' as const,
                genomicSourceClass: extractGenomicSourceClass(resource),
                variantAlleleFrequency: extractVariantAlleleFrequency(resource),
                dxImplications: [],
                txImplications: [],
                molecularConsequences: []
            } satisfies Variant)),
            ...extractVariantResources(structuralData).map(resource => {
                const structural = extractStructuralDetails(resource);
                return {
                    range,
                    resolvedRange: genomicRange,
                    sourceObservationId: resource.id,
                    variant: getStructuralDescription(structural),
                    variantType: 'structural' as const,
                    structural,
                    genomicSourceClass: extractGenomicSourceClass(resource),
                    variantAlleleFrequency: extractVariantAlleleFrequency(resource),
                    dxImplications: [],
                    txImplications: [],
                    molecularConsequences: []
                } satisfies Variant;
            })
        ];

    } catch (error) {
        console.error(`Error fetching variants for ${range}:`, error);
        throw error;
    }
};
