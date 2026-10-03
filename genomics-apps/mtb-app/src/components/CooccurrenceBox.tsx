import { useState } from 'react';
import { ChevronRight } from 'lucide-react';
import { CooccurrenceProfile, describeCooccurrenceEvidence, sortEvidenceByLevel } from '@/services/cooccurrenceService';

interface CooccurrenceBoxProps {
    profiles: CooccurrenceProfile[];
    defaultOpen?: boolean;
}

export default function CooccurrenceBox({ profiles, defaultOpen = false }: CooccurrenceBoxProps) {
    const [isOpen, setIsOpen] = useState(defaultOpen);

    if (profiles.length === 0) {
        return null;
    }

    return (
        <div className="mb-3 rounded-lg border border-violet-300 bg-violet-50">
            <button
                type="button"
                aria-expanded={isOpen}
                onClick={() => setIsOpen((currentValue) => !currentValue)}
                className="flex w-full items-center gap-1.5 px-3 py-2 text-left text-xs font-bold uppercase tracking-wide text-violet-700 focus:outline-none focus:ring-2 focus:ring-violet-300"
            >
                <ChevronRight className={`h-4 w-4 flex-shrink-0 transition-transform ${isOpen ? 'rotate-90' : ''}`} aria-hidden="true" />
                <span>Potentially relevant co-occurring variants</span>
                <span className="font-semibold normal-case tracking-normal text-gray-500">({profiles.length})</span>
            </button>
            {isOpen && (
                <div className="space-y-3 px-3 pb-3">
                    {profiles.map((profile) => (
                        <div key={profile.id} className="border-t border-violet-200 pt-2 first:border-t-0 first:pt-0">
                            <a
                                href={profile.url}
                                target="_blank"
                                rel="noopener noreferrer"
                                title="Open in CIViC"
                                className="text-sm font-semibold text-blue-600 hover:text-blue-800 hover:underline"
                            >
                                {profile.name}
                            </a>
                            <ul className="mt-1 list-disc space-y-1 pl-5 text-sm font-medium text-gray-900 marker:text-violet-500">
                                {sortEvidenceByLevel(profile.evidence).map((evidence) => (
                                    <li key={evidence.evidenceId}>{describeCooccurrenceEvidence(evidence)}</li>
                                ))}
                            </ul>
                        </div>
                    ))}
                </div>
            )}
        </div>
    );
}
