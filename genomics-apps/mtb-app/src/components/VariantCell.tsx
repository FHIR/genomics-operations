

interface VariantCellProps {
    variant?: string;
}

export default function VariantCell({ variant }: VariantCellProps) {
    if (!variant) {
        return (
            <td className="p-3">
                <span className="italic text-gray-400">Loading...</span>
            </td>
        );
    }

    // Split the variant string into genomic variant and protein change
    const parts = variant.split(' ');
    const genomicVariant = parts[0];
    const proteinChange = parts[1] ? parts[1].replace(/[()]/g, '') : null;

    return (
        <td className="p-3">
            <div className="text-gray-700">
                <span className="font-medium text-gray-900">{genomicVariant}</span>
                {proteinChange && (<>
                    <span className="mx-2 text-gray-500">(</span>
                    <span className="font-medium text-gray-900">{proteinChange}</span>
                    <span className="text-gray-500">)</span>
                </>)}
            </div>
        </td>
    );
}