import { useEffect, useState } from "react";
import { loadGeneListKbRows } from "@/cancerFilter/geneListKb";

interface CancerSelectProps {
    onSelect: (type: string) => void;
}

const CancerSelect = ({ onSelect }: CancerSelectProps) => {
    const [selectedType, setSelectedType] = useState<string>("");
    const [cancerTypes, setCancerTypes] = useState<string[]>([]);

    useEffect(() => {
        loadGeneListKbRows().then((rows) => {
            const availableTypes: string[] = Array.from(
                new Set(
                    rows
                        .map((row) => row["Cancer category"])
                        .filter((type): type is string => Boolean(type))
                )
            ).sort((left, right) => left.localeCompare(right));

            setCancerTypes(availableTypes);
        });
    }, []);

    const handleChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        const type = e.target.value;
        setSelectedType(type);
        onSelect(type);
    };

    return (
        <div>
            <label className="block mb-2 text-gray-600">
                Select Cancer Type
            </label>
            <select
                value={selectedType}
                onChange={handleChange}
                disabled={cancerTypes.length === 0}
                className="w-full p-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
                <option value="">{cancerTypes.length === 0 ? "Loading cancer types..." : "No cancer type selected"}</option>
                {cancerTypes.map((type) => (
                    <option key={type} value={type}>
                        {type}
                    </option>
                ))}
            </select>
        </div>
    );
};

export default CancerSelect;
