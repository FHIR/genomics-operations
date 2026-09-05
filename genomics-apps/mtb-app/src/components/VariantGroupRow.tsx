import {JSX, useState} from 'react';
import { Variant } from '@/types/variants';
import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import LoadingCell from './LoadingCell';

interface VariantGroupRowProps {
    range: string;
    variants: Variant[];
    renderDxImplications: (dx?: DxImplication[]) => JSX.Element;
    renderTxImplications: (tx?: ProcessedTxImplication[]) => JSX.Element;
    renderMolecularConsequences: (mc?: MolecularConsequence[]) => JSX.Element;
}

export default function VariantGroupRow({
                                            range,
                                            variants,
                                            renderDxImplications,
                                            renderTxImplications,
                                            renderMolecularConsequences,
                                        }: VariantGroupRowProps) {
    const [expanded, setExpanded] = useState(false);

    const visibleVariants = variants.length > 2 ? variants.slice(0, 2) : [variants[0]];
    const rest = variants.length > 2 ? variants.slice(2) : variants.slice(1);

    return (
        <>
            {visibleVariants.map((variant, index) => (
                <tr key={variant.id || `${range}-visible-${index}`} className="border-t border-gray-200 hover:bg-gray-50">
                    <td className="p-3">{index === 0 ? range : ''}</td>
                    <td className="p-3 max-w-0">
                        <div className="flex items-center gap-2">
                            <span className="truncate block">{variant.variant}</span>
                            {variant.isLoading && (
                                <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-500 flex-shrink-0"></div>
                            )}
                        </div>
                        {variant.molecularConsequences?.[0]?.proteinChange && (
                            <div className="text-sm text-gray-600 truncate">
                                ({variant.molecularConsequences[0].proteinChange})
                            </div>
                        )}
                    </td>
                    <td className="p-0">
                        <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                            {renderMolecularConsequences(variant.molecularConsequences)}
                        </LoadingCell>
                    </td>
                    <td className="p-0">
                        <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                            {renderDxImplications(variant.dxImplications)}
                        </LoadingCell>
                    </td>
                    <td className="p-0">
                        <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                            {renderTxImplications(variant.txImplications)}
                        </LoadingCell>
                    </td>
                </tr>
            ))}
            {expanded && rest.map((variant, i) => (
                <tr key={variant.id || `${range}-variant-${i}`} className="border-t border-gray-100 hover:bg-gray-50">
                    <td className="p-3"/>
                    <td className="p-3 max-w-0">
                        <div className="flex items-center gap-2">
                            <span className="truncate block">{variant.variant}</span>
                            {variant.isLoading && (
                                <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-500 flex-shrink-0"></div>
                            )}
                        </div>
                        {variant.molecularConsequences?.[0]?.proteinChange && (
                            <div className="text-sm text-gray-600 truncate">
                                ({variant.molecularConsequences[0].proteinChange})
                            </div>
                        )}
                    </td>
                    <td className="p-0">
                        <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                            {renderMolecularConsequences(variant.molecularConsequences)}
                        </LoadingCell>
                    </td>
                    <td className="p-0">
                        <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                            {renderDxImplications(variant.dxImplications)}
                        </LoadingCell>
                    </td>
                    <td className="p-0">
                        <LoadingCell isLoading={variant.isLoading} error={variant.error}>
                            {renderTxImplications(variant.txImplications)}
                        </LoadingCell>
                    </td>
                </tr>
            ))}
            {rest.length > 0 && (
                <tr>
                    <td colSpan={5} className="text-center p-2 text-sm">
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
