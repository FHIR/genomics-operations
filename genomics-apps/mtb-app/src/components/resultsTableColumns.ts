export type ColumnId =
    | 'range'
    | 'variant'
    | 'genomicSourceClass'
    | 'oncogenicityPrediction'
    | 'molecularConsequences'
    | 'dxImplications'
    | 'txImplications';

export interface ResultsTableColumnDefinition {
    id: ColumnId;
    label: string;
    defaultWidth: number;
    minWidth: number;
}

export const RESULTS_TABLE_COLUMNS: ResultsTableColumnDefinition[] = [
    { id: 'range', label: 'Range', defaultWidth: 220, minWidth: 180 },
    { id: 'variant', label: 'Variant', defaultWidth: 220, minWidth: 180 },
    { id: 'genomicSourceClass', label: 'Genomic Source Class', defaultWidth: 220, minWidth: 180 },
    {
        id: 'oncogenicityPrediction',
        label: 'Oncogenicity Prediction',
        defaultWidth: 220,
        minWidth: 180,
    },
    {
        id: 'molecularConsequences',
        label: 'Molecular Consequences',
        defaultWidth: 220,
        minWidth: 180,
    },
    { id: 'dxImplications', label: 'Dx Implications', defaultWidth: 220, minWidth: 180 },
    { id: 'txImplications', label: 'Tx Implications', defaultWidth: 150, minWidth: 150 },
];

export const DEFAULT_COLUMN_ORDER = RESULTS_TABLE_COLUMNS.map((column) => column.id);
