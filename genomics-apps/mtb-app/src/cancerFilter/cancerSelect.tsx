import { useState } from "react";

interface CancerSelectProps {
    onSelect: (type: string) => void;
}

const CancerSelect = ({ onSelect }: CancerSelectProps) => {
    const [selectedType, setSelectedType] = useState<string>("");

    const cancerTypes = [
        "Breast carcinoma",
        "NSCLC",
        // add more as needed
    ];

    const handleChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        const type = e.target.value;
        setSelectedType(type);
        onSelect(type);
    };

    return (
        <div className="mb-8">
            <label className="block mb-2 text-gray-600">
                Select Cancer Type
            </label>
            <select
                value={selectedType}
                onChange={handleChange}
                className="w-full p-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
                <option value="" disabled>Select a cancer type</option>
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
