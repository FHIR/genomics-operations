import {useState} from "react";
import dynamic from 'next/dynamic';

const TooltipIcon = dynamic(() => import('../components/TooltipIcon'), {
    ssr: false,
    loading: () => <span className="w-4 h-4 inline-block ml-1" />
});

interface ActionableCheckBoxesProps {
    onLabelChange: (label: string) => void;
}

const ACTIONABILITY_OPTIONS = [
    {
        id: 'actionable-this-tumor',
        label: 'Actionable, this tumor type',
        tooltip: 'Strong evidence for this cancer type'
    },
    {
        id: 'actionable-other-tumor',
        label: 'Actionable, other tumor type',
        tooltip: 'Strong evidence but for other cancer types'
    },
    {
        id: 'possibly-actionable',
        label: 'Possibly actionable',
        tooltip: 'Weaker or unvalidated evidence'
    },
    {
        id: 'actionable-none',
        label: 'None',
        tooltip: 'Do not apply an actionability filter'
    }
] as const;

export default function ActionableCheckBoxes({ onLabelChange }: ActionableCheckBoxesProps) {
    const [selectedLabel, setSelectedLabel] = useState<string>('None');

    const handleSelectionChange = (label: string) => {
        setSelectedLabel(label);
        onLabelChange(label === 'None' ? '' : label);
    };

    return (
        <div className="mt-4">
            <p className="text-gray-600 mb-2">Filter by actionability:</p>
            <div className="flex gap-4 lg:gap-8 xl:gap-12 2xl:gap-16 items-center flex-wrap">
                {ACTIONABILITY_OPTIONS.map((option) => (
                    <div key={option.id} className="flex items-center gap-2">
                        <input
                            id={option.id}
                            type="radio"
                            name="actionability-filter"
                            checked={selectedLabel === option.label}
                            onChange={() => handleSelectionChange(option.label)}
                            className="h-4 w-4 text-blue-500 border border-gray-400 focus:ring-2 focus:ring-blue-500 focus:ring-offset-0"
                        />
                        <label htmlFor={option.id} className="text-gray-600">
                            {option.label}
                        </label>
                        <TooltipIcon content={option.tooltip} />
                    </div>
                ))}
            </div>
        </div>
    );
}