export type ColumnId =
    | 'range'
    | 'variant'
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
    { id: 'range', label: 'Range', defaultWidth: 200, minWidth: 140 },
    { id: 'variant', label: 'Variant', defaultWidth: 410, minWidth: 220 },
    {
        id: 'oncogenicityPrediction',
        label: 'Oncogenicity Prediction',
        defaultWidth: 240,
        minWidth: 200,
    },
    {
        id: 'molecularConsequences',
        label: 'Molecular Consequences',
        defaultWidth: 280,
        minWidth: 220,
    },
    { id: 'dxImplications', label: 'Dx Implications', defaultWidth: 235, minWidth: 200 },
    { id: 'txImplications', label: 'Tx Implications', defaultWidth: 240, minWidth: 220 },
];

export const DEFAULT_COLUMN_ORDER = RESULTS_TABLE_COLUMNS.map((column) => column.id);
