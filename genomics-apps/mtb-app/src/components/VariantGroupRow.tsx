import { JSX, useEffect, useState } from 'react';
import { Variant } from '@/types/variants';
import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import { OncogenicityPredictionResult } from '@/types/oncogenicity';
import LoadingCell from './LoadingCell';
import OncogenicityPredictionCell from './OncogenicityPredictionCell';
import { ResultsTableColumnDefinition } from './resultsTableColumns';

interface VariantGroupRowProps {
    range: string;
    variants: Variant[];
    columns: ResultsTableColumnDefinition[];
    renderDxImplications: (dx?: DxImplication[]) => JSX.Element;
    renderTxImplications: (tx: ProcessedTxImplication[] | undefined, variant: Variant) => JSX.Element;
    renderMolecularConsequences: (mc?: MolecularConsequence[]) => JSX.Element;
    oncogenicityResults: Record<string, OncogenicityPredictionResult>;
    getOncogenicityKey: (variant: Variant) => string;
    onComputeOncogenicity: (variant: Variant) => void;
    onOpenOncogenicityDetails: (variant: Variant) => void;
    selectedVariantIds: Set<string>;
    onToggleVariantSelection: (variantId: string) => void;
    highlightedVariantId?: string | null;
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
    oncogenicityResults,
    getOncogenicityKey,
    onComputeOncogenicity,
    onOpenOncogenicityDetails,
    selectedVariantIds,
    onToggleVariantSelection,
    highlightedVariantId,
}: VariantGroupRowProps) {
    const [expanded, setExpanded] = useState(false);

    const visibleVariants = variants.length > 2 ? variants.slice(0, 2) : [variants[0]];
    const rest = variants.length > 2 ? variants.slice(2) : variants.slice(1);

    // Expand the group when a highlighted variant is in its collapsed rows
    useEffect(() => {
        if (highlightedVariantId && rest.some((variant) => variant.id === highlightedVariantId)) {
            setExpanded(true);
        }
    }, [highlightedVariantId, rest]);

    const getRowClassName = (variant: Variant, borderClassName: string) => {
        const isSelected = Boolean(variant.id && selectedVariantIds.has(variant.id));
        const isHighlighted = Boolean(variant.id && variant.id === highlightedVariantId);
        const backgroundClassName = isHighlighted ? 'bg-amber-100' : isSelected ? 'bg-blue-50/60 hover:bg-blue-50' : 'hover:bg-gray-50';

        return `border-t ${borderClassName} transition-colors duration-700 ${backgroundClassName}`;
    };

    const renderCells = (variant: Variant, showRange: boolean) => {
        const selectionCell = (
            <td key="select" className="p-3 align-top">
                {variant.id && (
                    <input
                        type="checkbox"
                        checked={selectedVariantIds.has(variant.id)}
                        onChange={() => onToggleVariantSelection(variant.id as string)}
                        aria-label={`Select ${range} ${variant.variant}`}
                        className="h-4 w-4 cursor-pointer accent-blue-600"
                    />
                )}
            </td>
        );

        return [selectionCell, ...columns.map((column) => {
            switch (column.id) {
                case 'range':
                    return <td key={column.id} className="p-3 align-top whitespace-normal break-all">{showRange ? range : ''}</td>;
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
                case 'genomicSourceClass':
                    return (
                        <td key={column.id} className="p-3 max-w-0 align-top">
                            <div className="block min-w-0 whitespace-normal break-words text-gray-700">
                                {variant.genomicSourceClass ?? '<unknown>'}
                            </div>
                        </td>
                    );
                case 'variantAlleleFrequency':
                    return (
                        <td key={column.id} className="p-3 max-w-0 align-top">
                            <div className="block min-w-0 whitespace-normal break-words text-gray-700">
                                {variant.variantAlleleFrequency ?? '<none found>'}
                            </div>
                        </td>
                    );
                case 'oncogenicityPrediction':
                    return (
                        <td key={column.id} className="p-3 align-top text-sm text-gray-600 whitespace-normal break-words">
                            <OncogenicityPredictionCell
                                variant={variant}
                                result={oncogenicityResults[getOncogenicityKey(variant)]}
                                onComputePrediction={onComputeOncogenicity}
                                onOpenDetails={onOpenOncogenicityDetails}
                            />
                        </td>
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
                                {renderTxImplications(variant.txImplications, variant)}
                            </LoadingCell>
                        </td>
                    );
            }
        })];
    };

    return (
        <>
            {visibleVariants.map((variant, index) => (
                <tr key={variant.id || `${range}-visible-${index}`} id={variant.id ? `variant-row-${variant.id}` : undefined} className={getRowClassName(variant, 'border-gray-200')}>
                    {renderCells(variant, index === 0)}
                </tr>
            ))}
            {expanded && rest.map((variant, i) => (
                <tr key={variant.id || `${range}-variant-${i}`} id={variant.id ? `variant-row-${variant.id}` : undefined} className={getRowClassName(variant, 'border-gray-100')}>
                    {renderCells(variant, false)}
                </tr>
            ))}
            {rest.length > 0 && (
                <tr>
                    <td colSpan={columns.length + 1} className="text-center p-2 text-sm">
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
