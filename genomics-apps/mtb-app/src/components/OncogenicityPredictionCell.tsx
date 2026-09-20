import Image from 'next/image';
import { Variant } from '@/types/variants';
import { getGaugeImagePath, OncogenicityPredictionResult } from '@/types/oncogenicity';

interface OncogenicityPredictionCellProps {
    variant: Variant;
    result?: OncogenicityPredictionResult;
    onComputePrediction: (variant: Variant) => void;
    onOpenDetails: (variant: Variant) => void;
}

export default function OncogenicityPredictionCell({
    variant,
    result,
    onComputePrediction,
    onOpenDetails,
}: OncogenicityPredictionCellProps) {
    const hasResolvedScore = typeof result?.score === 'number' && !Number.isNaN(result.score);

    if (variant.variantType !== 'simple') {
        return (
            <div className="text-xs leading-5 text-gray-500 whitespace-normal break-words">
                Currently only applicable to simple variants
            </div>
        );
    }

    if (!result || (result.status === 'ready' && !hasResolvedScore && result.gauge === 'undetermined')) {
        return (
            <button
                type="button"
                onClick={() => onComputePrediction(variant)}
                className="rounded-md border border-blue-700 bg-blue-600 px-3 py-2 text-xs font-semibold text-white shadow-[0_2px_0_0_rgb(29_78_216)] transition-[transform,box-shadow,background-color] hover:bg-blue-700 hover:shadow-[0_1px_0_0_rgb(30_64_175)] active:translate-y-px active:shadow-none focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
            >
                Compute prediction
            </button>
        );
    }

    if (result.status === 'loading') {
        return (
            <div className="flex items-center gap-2 text-xs font-medium text-blue-700">
                <span className="h-4 w-4 animate-spin rounded-full border-2 border-blue-200 border-t-blue-600" />
                <span>Computing...</span>
            </div>
        );
    }

    return (
        <button
            type="button"
            onClick={() => onOpenDetails(variant)}
            className="flex flex-col items-start gap-2 rounded-lg bg-white p-2 shadow-sm ring-1 ring-slate-200 transition-shadow hover:shadow-md focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
        >
            <Image
                src={getGaugeImagePath(result.gauge)}
                alt="Oncogenicity prediction gauge"
                width={128}
                height={72}
                className="h-auto w-32"
            />
        </button>
    );
}
