import { DxImplication } from '@/services/dxService';

interface DxImplicationCellProps {
    implications?: DxImplication[];
}

export default function DxImplicationCell({ implications }: DxImplicationCellProps) {
    if (!implications || implications.length === 0) {
        return (
            <div className="p-3 text-gray-400 italic">&lt;none found&gt;</div>
        );
    }

    function evidenceLevelToStars(level: string): number {
        switch (level.toLowerCase()) {
            case 'practice guideline':
                return 4;
            case 'reviewed by expert panel':
                return 3;
            case 'criteria provided, single submitter':
                return 1;
            case 'no assertion criteria provided':
            case 'no classification provided':
                return 0;
            default:
                return 0;
        }
    }

    function pickBestDxImplication(implications: DxImplication[]): DxImplication | null {
        const concordanceGroups = [
            ['Benign', 'Likely benign'],
            ['Pathogenic', 'Likely pathogenic', 'Likely pathogenic, low penetrance', 'Pathogenic, low penetrance'],
        ];

        const withStars = implications.map((imp) => ({
            imp,
            stars: imp.evidenceLevel ? evidenceLevelToStars(imp.evidenceLevel) : 0,
        }));

        const fourStar = withStars.find(item => item.stars === 4);
        if (fourStar) return fourStar.imp;

        const threeStar = withStars.find(item => item.stars === 3);
        if (threeStar) return threeStar.imp;

        const oneStars = withStars.filter(item => item.stars === 1);
        if (oneStars.length > 1) {
            const sigs = oneStars.map(item => item.imp.clinicalSignificance?.toLowerCase() || '');
            const isConcordant = concordanceGroups.some(group =>
                sigs.every(sig => group.map(g => g.toLowerCase()).includes(sig))
            );
            // If implications are concordant, return the first one, otherwise null
            return isConcordant ? oneStars[0].imp : null;
        }
        if (oneStars.length === 1) {
            return oneStars[0].imp;
        }

        const zeroStars = withStars.filter(item => item.stars === 0);
        if (zeroStars.length > 0) {
            return zeroStars[0].imp;
        }

        return null;
    }

    const bestImplication = pickBestDxImplication(implications);

    if (!bestImplication) {
        return (
            <div className="p-3 text-gray-400 italic">&lt;none found&gt;</div>
        );
    }

    const stars = evidenceLevelToStars(bestImplication.evidenceLevel || '');

    return (
        <div className="p-3 text-gray-700 font-medium">
            {bestImplication.clinicalSignificance || 'Unknown'}{' '}
            <a
                href={bestImplication.clinvarLink}
                target="_blank"
                rel="noopener noreferrer"
                className="text-blue-600 hover:underline"
            >
                (ClinVar; {stars} {stars === 1 ? 'star' : 'stars'})
            </a>
        </div>
    );
}
