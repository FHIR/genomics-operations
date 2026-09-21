import { useState } from 'react';

interface SearchFormProps {
  searchInput: string;
  setSearchInput: (value: string) => void;
  handleSearch: () => void;
  selectedCancerType: string;
  onPresetSelect: (label: string) => void;
  enableCatVrsQueries: boolean;
  onEnableCatVrsQueriesChange: (enabled: boolean) => void;
  embedded?: boolean;
  className?: string;
}

interface MrnSelectorProps {
  subjectId: string;
  setSubjectId: (value: string) => void;
  className?: string;
}

const patientIDs: string[] = [
  "ABC123", "NA19240", "ABC456",
  "NA18498", "NA19247", "m123",
  "ABC789", "NA18499", "NA19256",
  "CA12345", "NA18870", "NB6TK328",
  "HCC1143", "NA18871", "NB6TK329",
  "HG00403", "NA19190", "XYZ123",
  "HG00406", "NA19210", "XYZ234",
  "HG02657", "NA19238", "XYZ345",
  "huC30902", "NA19239", "L2345"
].sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));

export function MrnSelector({
  subjectId,
  setSubjectId,
  className = 'mb-8',
}: MrnSelectorProps) {
  const [isCustomMRN, setIsCustomMRN] = useState(false);
  const [customMRN, setCustomMRN] = useState('');

  const handleMRNChange = (value: string) => {
    if (value === 'custom') {
      setIsCustomMRN(true);
      setSubjectId('');
    } else {
      setIsCustomMRN(false);
      setSubjectId(value);
      setCustomMRN('');
    }
  };

  const handleCustomMRNChange = (value: string) => {
    setCustomMRN(value);
    setSubjectId(value);
  };

  return (
    <div className={className.trim()}>
      <div>
        <label htmlFor="mrn-select" className="block mb-2 text-gray-600">MRN (Medical Record Number)</label>
        {!isCustomMRN ? (
          <div className="flex gap-2">
            <select
              id="mrn-select"
              value={subjectId}
              onChange={(e) => handleMRNChange(e.target.value)}
              className="w-full p-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
              <option value="">Select a patient MRN</option>
              {patientIDs.map((id) => (
                <option key={id} value={id}>
                  {id}
                </option>
              ))}
              <option value="custom">Enter custom MRN...</option>
            </select>
          </div>
        ) : (
          <div className="flex gap-2">
            <input
              type="text"
              value={customMRN}
              onChange={(e) => handleCustomMRNChange(e.target.value)}
              placeholder="Enter custom MRN"
              className="w-full p-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            />
            <button
              onClick={() => {
                setIsCustomMRN(false);
                setSubjectId('');
                setCustomMRN('');
              }}
              className="rounded-md border border-gray-300 px-4 py-2 text-sm text-gray-700 transition-colors hover:bg-gray-50"
            >
              Back to List
            </button>
          </div>
        )}
      </div>
    </div>
  );
}

export default function SearchForm({
  searchInput,
  setSearchInput,
  handleSearch,
  selectedCancerType,
  onPresetSelect,
  enableCatVrsQueries,
  onEnableCatVrsQueriesChange,
  embedded = false,
  className = 'mb-8',
}: SearchFormProps) {
  const presetButtons = [
    {
      label: 'Actionable Genes',
      value: 'Actionable Gene List',
    },
    {
      label: 'Extended Gene List',
      value: 'Extended Gene List',
    },
  ] as const;

  const containerClassName = embedded
    ? className.trim()
    : `bg-gray-100 p-8 rounded-lg shadow-sm ${className}`.trim();

  const presetHelpText = selectedCancerType
    ? 'Available for the selected cancer type.'
    : 'Select a cancer type to enable.';

  return (
    <div className={containerClassName}>
      <div className="rounded-xl border border-blue-100 bg-white p-4 shadow-[inset_0_1px_0_0_rgba(255,255,255,0.9)]">
        <div className="mb-3 flex items-center justify-between gap-3">
          <div>
            <label htmlFor="search-terms-input" className="inline-flex items-center gap-2 text-sm font-semibold text-slate-900">
              <span className="inline-block h-2.5 w-2.5 rounded-full bg-sky-500 shadow-[0_0_0_4px_rgba(14,165,233,0.14)]" />
              <span>Search terms</span>
            </label>
            <p className="mt-1 text-sm text-slate-600">Gene symbols or genomic ranges, separated by commas. Ranges must use zero-based RefSeq:start-end format.</p>
          </div>
          <button
            onClick={handleSearch}
            className="shrink-0 rounded-md border border-blue-700 bg-blue-600 px-5 py-2.5 text-sm font-semibold text-white shadow-[0_3px_0_0_rgb(29_78_216)] transition-[transform,box-shadow,background-color] hover:bg-blue-700 hover:shadow-[0_2px_0_0_rgb(30_64_175)] active:translate-y-px active:shadow-[0_1px_0_0_rgb(30_64_175)] focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
          >
            Search
          </button>
        </div>

        <input
          id="search-terms-input"
          type="text"
          value={searchInput}
          onChange={(e) => setSearchInput(e.target.value)}
          className="w-full rounded-xl border-2 border-sky-300 bg-sky-50/55 px-4 py-4 text-lg text-slate-950 shadow-[inset_0_1px_0_0_rgba(255,255,255,0.95),0_2px_10px_rgba(14,116,144,0.08)] outline-none transition-[border-color,box-shadow,background-color] placeholder:text-slate-500 focus:border-sky-500 focus:bg-white focus:ring-4 focus:ring-sky-100"
          onKeyDown={(e) => e.key === 'Enter' && handleSearch()}
          placeholder="BRAF, EGFR, ALK"
        />

        <div className="mt-3 border-t border-gray-100 pt-3">
          <div className="flex flex-wrap items-center gap-2 text-xs text-gray-500">
            <span className="font-semibold text-gray-700">Examples:</span>
            <span className="font-mono">BRAF</span>
            <span className="text-gray-300">•</span>
            <span className="font-mono">NC_000007.14:140713327-140924929</span>
            <span className="text-gray-300">•</span>
            <span className="font-mono">NC_000007.14:55174721-55174820, BRAF, EGFR, ALK</span>
          </div>

          <div className="mt-3 flex flex-wrap items-center gap-2 text-xs text-gray-500">
            <span className="font-semibold text-gray-700">Presets:</span>
            <span>{presetHelpText}</span>
            {presetButtons.map((preset) => {
              return (
                <button
                  key={preset.value}
                  type="button"
                  onClick={() => onPresetSelect(preset.value)}
                  disabled={!selectedCancerType}
                  className="inline-flex items-center rounded-full border border-slate-300 bg-slate-50 px-3 py-1.5 text-xs font-semibold text-slate-700 transition-colors hover:border-blue-300 hover:bg-blue-50 hover:text-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2 disabled:cursor-not-allowed disabled:border-gray-300 disabled:bg-gray-100 disabled:text-gray-400"
                >
                  {preset.label}
                </button>
              );
            })}
          </div>
        </div>
      </div>

      <div className="mt-3 flex justify-end">
        <label className="flex items-center gap-2 rounded-md border border-gray-200 bg-white px-3 py-2 text-sm text-gray-700 shadow-sm">
          <input
            type="checkbox"
            checked={enableCatVrsQueries}
            onChange={(event) => onEnableCatVrsQueriesChange(event.target.checked)}
          />
          <span className="font-medium">Enable Cat-VRS queries</span>
        </label>
      </div>
    </div>
  );
}
