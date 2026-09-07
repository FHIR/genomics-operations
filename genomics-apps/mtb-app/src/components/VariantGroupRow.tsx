import { JSX, useState } from 'react';
import { Variant } from '@/types/variants';
import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import LoadingCell from './LoadingCell';
import { ResultsTableColumnDefinition } from './resultsTableColumns';

interface VariantGroupRowProps {
    range: string;
    variants: Variant[];
    columns: ResultsTableColumnDefinition[];
    renderDxImplications: (dx?: DxImplication[]) => JSX.Element;
    renderTxImplications: (tx?: ProcessedTxImplication[]) => JSX.Element;
    renderMolecularConsequences: (mc?: MolecularConsequence[]) => JSX.Element;
}

const getSimpleVariantLabel = (variantString: string) => {
    const parts = variantString.split(':');

    if (parts.length !== 4) {
        return 'Simple';
    }

    const [, , deleted = '', inserted = ''] = parts;

    if (deleted.length === inserted.length) {
        if (deleted.length === 1) {
            return 'SNV';
        }

        if (deleted.length > 1) {
            return 'MNV';
        }
    }

    if (deleted.length !== inserted.length) {
        return 'InDel';
    }

    return 'Simple';
};

const renderVariantLabel = (variant: Variant) => {
    if (variant.variantType === 'simple') {
        const simpleVariantLabel = getSimpleVariantLabel(variant.variant);

        return (
            <>
                <span className="font-semibold">{simpleVariantLabel}</span>{' '}
                <span>({variant.variant})</span>
            </>
        );
    }

    const structuralMatch = /^(.*?)(\s\([^)]+\))?(\s\(Copies:.*\))?$/.exec(variant.variant);

    if (!structuralMatch) {
        return <strong>{variant.variant}</strong>;
    }

    const [, dnaChangeType = variant.variant, locationSuffix = '', copiesSuffix = ''] = structuralMatch;

    return (
        <>
            <span className="font-semibold">{dnaChangeType}</span>
            <span>{locationSuffix}{copiesSuffix}</span>
        </>
    );
};

export default function VariantGroupRow({
    range,
    variants,
    columns,
    renderDxImplications,
    renderTxImplications,
    renderMolecularConsequences,
}: VariantGroupRowProps) {
    const [expanded, setExpanded] = useState(false);

    const visibleVariants = variants.length > 2 ? variants.slice(0, 2) : [variants[0]];
    const rest = variants.length > 2 ? variants.slice(2) : variants.slice(1);

    const renderCells = (variant: Variant, showRange: boolean) => {
        return columns.map((column) => {
            switch (column.id) {
                case 'range':
                    return <td key={column.id} className="p-3">{showRange ? range : ''}</td>;
                case 'variant':
                    return (
                        <td key={column.id} className="p-3 max-w-0 align-top">
                            <div className="flex min-w-0 items-start gap-2">
                                <span className="block min-w-0 whitespace-normal break-words">{renderVariantLabel(variant)}</span>
                                {variant.isLoading && (
                                    <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-500 flex-shrink-0"></div>
                                )}
                            </div>
                            {variant.molecularConsequences?.[0]?.proteinChange && (
                                <div className="mt-1 text-sm text-gray-600 whitespace-normal break-words">
                                    ({variant.molecularConsequences[0].proteinChange})
                                </div>
                            )}
                        </td>
                    );
                case 'oncogenicityPrediction':
                    return (
                        <td key={column.id} className="p-3 align-top text-sm text-gray-600" />
                    );
                case 'molecularConsequences':
                    return (
                        <td key={column.id} className="p-0 align-top">
                            <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                                {renderMolecularConsequences(variant.molecularConsequences)}
                            </LoadingCell>
                        </td>
                    );
                case 'dxImplications':
                    return (
                        <td key={column.id} className="p-0 align-top">
                            <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                                {renderDxImplications(variant.dxImplications)}
                            </LoadingCell>
                        </td>
                    );
                case 'txImplications':
                    return (
                        <td key={column.id} className="p-0 align-top">
                            <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                                {renderTxImplications(variant.txImplications)}
                            </LoadingCell>
                        </td>
                    );
            }
        });
    };

    return (
        <>
            {visibleVariants.map((variant, index) => (
                <tr key={variant.id || `${range}-visible-${index}`} className="border-t border-gray-200 hover:bg-gray-50">
                    {renderCells(variant, index === 0)}
                </tr>
            ))}
            {expanded && rest.map((variant, i) => (
                <tr key={variant.id || `${range}-variant-${i}`} className="border-t border-gray-100 hover:bg-gray-50">
                    {renderCells(variant, false)}
                </tr>
            ))}
            {rest.length > 0 && (
                <tr>
                    <td colSpan={columns.length} className="text-center p-2 text-sm">
                        <button
                            className="text-blue-600 hover:underline"
                            onClick={() => setExpanded(!expanded)}
                        >
                            {expanded ? 'Hide additional variants' : `Show ${rest.length} more variant(s)`}
                        </button>
                    </td>
                </tr>
            )}
        </>
    );
}
