import { useState } from 'react';
import { ProcessedTxImplication } from '@/services/txService';
import { useMemo } from 'react';

interface TxImplicationCellProps {
    implications?: ProcessedTxImplication[];
}

function ImplicationDisplay({ implication }: { implication: ProcessedTxImplication }) {
    // Build the display string in the format: "Variant {Implication} to {Medication} in {txPhenotype} [{txEvidence}]"
    const buildDisplayString = () => {
        const parts = [];

        // Start with "Variant"
        parts.push("Variant");

        // Add therapeutic implication if available
        if (implication.therapeuticImplication) {
            parts.push(implication.therapeuticImplication);
        }

        // Add "to" and medication if available
        if (implication.medication) {
            parts.push("to", implication.medication);
        }

        // Add "in" and phenotypic context if available
        if (implication.phenotypicContext) {
            parts.push("in", implication.phenotypicContext);
        }

        // Add evidence level in brackets if available
        if (implication.evidenceLevel) {
            parts.push(`[${implication.evidenceLevel}]`);
        }

        return parts.join(" ");
    };

    return (
        <div className="mb-2 border-t border-gray-200 pt-2">
            <div className="text-gray-700 mb-2">
                <span className="font-medium text-gray-900">{buildDisplayString()}</span>
            </div>
            {implication.hyperlink && (
                <a
                    href={implication.hyperlink}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 hover:text-blue-800 hover:underline text-sm"
                >
                    {implication.hyperlink.includes('civicdb.org/variants/')
                        ? 'CIViC entry'
                        : implication.hyperlink.includes('ncbi.nlm.nih.gov/clinvar/variation/')
                        ? 'ClinVar entry'
                        : implication.hyperlink.includes('cancer.gov')
                        ? 'NCI PDQ entry'
                        : 'External source'}
                </a>
            )}
        </div>
    );
}

export default function TxImplicationCell({ implications }: TxImplicationCellProps) {
    const [viewMode, setViewMode] = useState<'limited' | 'distinct' | 'all'>('limited'); // Default to showing limited (first 2 distinct)

    // Calculate distinct implications
    const distinctImplications = useMemo(() => {
        if (!implications) return [];
        const cache = new Set();
        return implications.filter((item) => {
            const key = [
                item.evidenceLevel,
                item.medication,
                item.phenotypicContext,
                item.therapeuticImplication
            ].join("|");
            if (cache.has(key)) return false;
            cache.add(key);
            return true;
        });
    }, [implications]);

    // Determine which implications to show based on view mode
    const implicationsToShow = useMemo(() => {
        switch (viewMode) {
            case 'limited':
                return distinctImplications.slice(0, 2);
            case 'distinct':
                return distinctImplications;
            case 'all':
                return implications || [];
            default:
                return distinctImplications.slice(0, 2);
        }
    }, [viewMode, distinctImplications, implications]);

    if (!implications) {
        return (
            <div className="p-3">
                <span className="italic text-gray-400">Loading...</span>
            </div>
        );
    }

    if (implications.length === 0) {
        return (
            <div className="p-3">
                <span className="text-gray-400">&lt;none found&gt;</span>
            </div>
        );
    }

    return (
        <div className="p-3">
            {/* Show implications based on current mode */}
            {implicationsToShow.map((implication, index) => (
                <ImplicationDisplay key={`${viewMode}-${index}`} implication={implication} />
            ))}

            {/* Show toggle buttons if there are more than 2 distinct implications or if distinct differs from all */}
            {(distinctImplications.length > 2 || distinctImplications.length !== implications.length) && (
                <div className="flex gap-4 mt-2">
                    {/* Show "View all distinct" button only if we're in limited mode and there are more than 2 distinct */}
                    {viewMode === 'limited' && distinctImplications.length > 2 && (
                        <button
                            onClick={() => setViewMode('distinct')}
                            className="text-sm text-blue-600 hover:underline"
                        >
                            View all distinct implications ({distinctImplications.length})
                        </button>
                    )}

                    {/* Show "View less" button if we're showing all distinct and there are more than 2 */}
                    {viewMode === 'distinct' && distinctImplications.length > 2 && (
                        <button
                            onClick={() => setViewMode('limited')}
                            className="text-sm text-blue-600 hover:underline"
                        >
                            View less
                        </button>
                    )}

                    {/* Show "View all implications" button if distinct differs from all and we're in limited mode */}
                    {distinctImplications.length !== implications.length && viewMode === 'limited' && (
                        <button
                            onClick={() => setViewMode('all')}
                            className="text-sm text-blue-600 hover:underline"
                        >
                            View all implications ({implications.length})
                        </button>
                    )}

                    {/* Show "View less" button if we're showing all implications */}
                    {viewMode === 'all' && (
                        <button
                            onClick={() => setViewMode('limited')}
                            className="text-sm text-blue-600 hover:underline"
                        >
                            View less
                        </button>
                    )}
                </div>
            )}
        </div>
    );
}