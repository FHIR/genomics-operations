import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';

/**
 * A structural variant's or fusion's fields, as recorded: the genomics server's rows and the reports give them
 * as separate components, not as HGVS. Every field is absent when not recorded.
 */
export interface StructuralDetails {
    changeType?: string; // The DNA change type's Sequence Ontology display, e.g. "copy_number_variation"
    referenceSequence?: string; // e.g. "NC_000006.11"
    start?: number; // Outer (else inner) range, 0-based
    end?: number;
    copyNumber?: number;
    insertedSequence?: string; // An insertion's literal inserted bases
    fusionGenes?: string[]; // A fusion's partner genes, in recorded order (5' partner first); absent for other types
}

// A genetic report a variant or implication came from (for the Reported chip and Report tags)
export interface ReportCredit {
    title: string;
    date?: string;
}

export interface Variant {
    id?: string;  // Unique identifier for the variant
    range: string;  // The original search term (gene or range)
    resolvedRange?: string; // The normalized genomic range used for API calls
    sourceObservationId?: string; // Back-reference to the source FHIR observation
    variant: string; // The variant identifier; for a structural variant or fusion, its description
    variantType?: 'simple' | 'structural';
    structural?: StructuralDetails; // A structural variant's or fusion's fields
    genomicSourceClass?: string;
    variantAlleleFrequency?: number; // Sample VAF (LOINC 81258-6)
    oncogenicityPrediction?: string; // Optional oncogenicity assessment when available
    dxImplications: DxImplication[]; // Diagnostic implications
    txImplications: ProcessedTxImplication[]; // Therapeutic implications
    molecularConsequences: MolecularConsequence[]; // Molecular consequences
    isLoading?: boolean; // Optional loading state for incremental display
    searchId?: number; // Optional search ID for tracking which search operation this belongs to
    error?: string; // Optional error message if fetching failed
    reportedIn?: ReportCredit[]; // The patient's genetic reports that include this variant
    isReportOnly?: boolean; // Built from a report; the genomics server doesn't have this variant
}
