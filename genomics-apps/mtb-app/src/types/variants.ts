import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';

export interface Variant {
    id?: string;  // Unique identifier for the variant
    range: string;  // The original search term (gene or range)
    resolvedRange?: string; // The normalized genomic range used for API calls
    sourceObservationId?: string; // Back-reference to the source FHIR observation
    variant: string; // The variant identifier
    variantType?: 'simple' | 'structural';
    genomicSourceClass?: string;
    oncogenicityPrediction?: string; // Optional oncogenicity assessment when available
    dxImplications: DxImplication[]; // Diagnostic implications
    txImplications: ProcessedTxImplication[]; // Therapeutic implications
    molecularConsequences: MolecularConsequence[]; // Molecular consequences
    isLoading?: boolean; // Optional loading state for incremental display
    searchId?: number; // Optional search ID for tracking which search operation this belongs to
    error?: string; // Optional error message if fetching failed
}
