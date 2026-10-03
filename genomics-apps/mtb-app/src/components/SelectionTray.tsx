import { Variant } from '@/types/variants';
import { getVariantGene } from '@/lib/variantDisplay';

interface SelectionTrayProps {
    selectedVariants: Variant[];
    hiddenCount: number;
    message?: string | null;
    onClear: () => void;
    onViewPathways: () => void;
}

export default function SelectionTray({ selectedVariants, hiddenCount, message, onClear, onViewPathways }: SelectionTrayProps) {
    if (selectedVariants.length === 0) {
        return null;
    }

    const geneCounts = selectedVariants.reduce<Record<string, number>>((counts, variant) => {
        const label = getVariantGene(variant) ?? variant.range;
        counts[label] = (counts[label] ?? 0) + 1;
        return counts;
    }, {});
    const genes = Object.entries(geneCounts);

    return (
        <div className="fixed inset-x-0 bottom-0 z-40 border-t-2 border-blue-600 bg-white px-6 py-3 shadow-[0_-6px_20px_rgba(0,0,0,0.12)]">
            <div className="mx-auto flex max-w-[1400px] flex-wrap items-center gap-x-4 gap-y-2">
                <span className="font-semibold text-gray-900">
                    {selectedVariants.length} variant{selectedVariants.length === 1 ? '' : 's'} · {genes.length} gene{genes.length === 1 ? '' : 's'}
                </span>
                <span className="flex min-w-0 flex-1 flex-wrap gap-1.5">
                    {genes.map(([gene, count]) => (
                        <span key={gene} className="rounded-md border border-gray-200 bg-gray-50 px-2 py-0.5 text-xs font-semibold text-gray-700">
                            {gene}{count > 1 && <span className="font-medium text-gray-500"> ×{count}</span>}
                        </span>
                    ))}
                </span>
                {hiddenCount > 0 && (
                    <span className="text-xs font-medium text-amber-700">{hiddenCount} hidden by current filters</span>
                )}
                {message && <span className="text-xs font-medium text-amber-700">{message}</span>}
                <button
                    type="button"
                    onClick={onClear}
                    className="rounded-md px-3 py-2 text-sm font-medium text-gray-700 hover:bg-gray-100 focus:outline-none focus:ring-2 focus:ring-blue-300"
                >
                    Clear
                </button>
                <button
                    type="button"
                    onClick={onViewPathways}
                    className="inline-flex items-center gap-2 rounded-md border border-blue-700 bg-blue-600 px-4 py-2.5 text-sm font-semibold text-white shadow-[0_3px_0_0_rgb(29_78_216)] hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
                >
                    View on pathways →
                </button>
            </div>
        </div>
    );
}
