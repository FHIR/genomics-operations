import { Variant } from '@/types/variants';

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
    component?: {
        code?: {
            coding?: {
                code: string;
            }[];
        };
        valueCodeableConcept?: {
            coding?: {
                code: string;
            }[];
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

        const response = await fetch(
            `${API_CONFIG.baseUrl}${API_CONFIG.endpoints.subjectOps}/$find-subject-variants?subject=${subjectId}&ranges=${encodeURIComponent(genomicRange)}&includeVariants=true`
        );

        if (!response.ok) {
            throw new Error(`API error: ${response.status} ${response.statusText}`);
        }

        const data = await response.json() as FhirResponse;

        // Process FHIR response to extract variants
        const variantResources = data.parameter
            ?.find(p => p.name === "variants")
            ?.part
            ?.filter(p => p.name === "variant" && p.resource)
            ?.map(p => p.resource) || [];

        // Create variant objects with required properties
        return variantResources.map(resource => {
            if (!resource) {
                throw new Error('Invalid variant resource');
            }

            return {
                range,
                variant: extractVariantString(resource),
                dxImplications: [], // Initialize with empty arrays
                txImplications: [],
                molecularConsequences: []
            } satisfies Variant;
        });

    } catch (error) {
        console.error(`Error fetching variants for ${range}:`, error);
        throw error;
    }
};