import { MolecularConsequence } from '@/services/mcService';

interface MolecularConsequenceCellProps {
    consequences?: MolecularConsequence[];
}

const severityOrder = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER'];

function getSeverityRank(impact?: string): number {
    return impact ? severityOrder.indexOf(impact.toUpperCase()) : severityOrder.length;
}

export default function MolecularConsequenceCell({ consequences }: MolecularConsequenceCellProps) {
    if (!consequences) {
        return (
            <div className="p-3">
                <span className="italic text-gray-400">Loading...</span>
            </div>
        );
    }

    if (consequences.length === 0) {
        return (
            <div className="p-3">
                <span className="text-gray-400">&lt;none found&gt;</span>
            </div>
        );
    }

    // Sort by impact severity and pick the first most severe consequence
    const sorted = [...consequences].sort((a, b) => getSeverityRank(a.impact) - getSeverityRank(b.impact));
    const mostSevere = sorted[0];

    return (
        <div className="p-3">
            {mostSevere.featureConsequence && (
                <div className="text-gray-700">
                    <span className="font-medium text-gray-900">{mostSevere.impact}</span>
                    <span className="mx-2 text-gray-500">/</span>
                    <span className="font-medium text-gray-900">{mostSevere.featureConsequence}</span>
                    {mostSevere.hyperlink && (
                        <span className="text-gray-500">
                            {' '}(
                            <a
                                href={mostSevere.hyperlink}
                                target="_blank"
                                rel="noopener noreferrer"
                                className="text-blue-600 hover:underline"
                            >
                                SnpEff
                            </a>
                            )
                        </span>
                    )}
                </div>
            )}
        </div>
    );
}