import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';

export interface Variant {
    id?: string;  // Unique identifier for the variant
    range: string;  // The original search term (gene or range)
    variant: string; // The variant identifier
    oncogenicityPrediction?: string; // Optional oncogenicity assessment when available
    dxImplications: DxImplication[]; // Diagnostic implications
    txImplications: ProcessedTxImplication[]; // Therapeutic implications
    molecularConsequences: MolecularConsequence[]; // Molecular consequences
    isLoading?: boolean; // Optional loading state for incremental display
    searchId?: number; // Optional search ID for tracking which search operation this belongs to
    error?: string; // Optional error message if fetching failed
}
