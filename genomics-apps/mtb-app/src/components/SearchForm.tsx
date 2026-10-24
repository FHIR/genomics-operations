import { useEffect, useState } from 'react';
import {
  getCancerSpecificPresetLabels,
  getGlobalGeneListLabels,
  GLOBAL_GENE_LIST_CATEGORY,
  loadGeneListKbRows,
  PresetSelection,
} from '@/cancerFilter/geneListKb';

interface SearchFormProps {
  searchInput: string;
  setSearchInput: (value: string) => void;
  handleSearch: () => void;
  // While a search runs, the button shows progress and can't start another
  isSearching: boolean;
  selectedCancerType: string;
  onPresetSelect: (selection: PresetSelection) => void;
  enableCatVrsQueries: boolean;
  onEnableCatVrsQueriesChange: (enabled: boolean) => void;
}

const presetClassName = 'rounded-full px-3 py-1 text-sm font-medium ring-1 ring-inset transition-colors focus:outline-none focus:ring-2 focus:ring-blue-300';

export default function SearchForm({
  searchInput,
  setSearchInput,
  handleSearch,
  isSearching,
  selectedCancerType,
  onPresetSelect,
  enableCatVrsQueries,
  onEnableCatVrsQueriesChange,
}: SearchFormProps) {
  const [cancerPresetLabels, setCancerPresetLabels] = useState<string[]>([]);
  const [globalGeneListLabels, setGlobalGeneListLabels] = useState<string[]>([]);

  useEffect(() => {
    let isMounted = true;

    loadGeneListKbRows().then((rows) => {
      if (!isMounted) {
        return;
      }

      setGlobalGeneListLabels(getGlobalGeneListLabels(rows));
      setCancerPresetLabels(getCancerSpecificPresetLabels(rows, selectedCancerType));
    });

    return () => {
      isMounted = false;
    };
  }, [selectedCancerType]);

  return (
    <div className="flex flex-col gap-3">
      <div className="flex flex-col gap-2.5 sm:flex-row">
        <label htmlFor="search-terms-input" className="sr-only">Search terms</label>
        <input
          id="search-terms-input"
          type="text"
          value={searchInput}
          onChange={(e) => setSearchInput(e.target.value)}
          onKeyDown={(e) => e.key === 'Enter' && !isSearching && handleSearch()}
          placeholder="BRAF, EGFR, ALK"
          className="min-w-0 flex-1 rounded-xl border-[1.5px] border-slate-400 bg-white px-4 py-3.5 text-base text-slate-900 shadow-sm outline-none transition-[border-color,box-shadow] placeholder:text-slate-400 focus:border-emerald-600 focus:ring-4 focus:ring-emerald-100"
        />
        <button
          type="button"
          onClick={handleSearch}
          disabled={isSearching}
          className="shrink-0 rounded-xl bg-slate-900 px-7 py-3 text-base font-semibold text-white transition-colors hover:bg-slate-800 focus:outline-none focus:ring-2 focus:ring-slate-400 focus:ring-offset-2 disabled:cursor-not-allowed disabled:bg-slate-400"
        >
          {isSearching ? 'Searching…' : 'Search'}
        </button>
      </div>

      <p className="-mt-1 text-sm text-slate-500">
        Gene symbols or genomic ranges, separated by commas. Ranges use zero-based RefSeq:start-end format, e.g.{' '}
        <code className="rounded bg-slate-100 px-1.5 py-0.5 font-mono text-xs text-slate-700">NC_000007.14:55174721-55174820</code>.
      </p>

      <div className="flex flex-wrap items-center gap-2 text-sm">
        <span className="mr-0.5 font-semibold text-slate-700">Presets</span>
        {/* Cancer-specific presets follow the Dx in the banner; general lists are always available */}
        {cancerPresetLabels.length > 0 ? (
          <>
            <span className="text-slate-500">{selectedCancerType}:</span>
            {cancerPresetLabels.map((label) => (
              <button
                key={`cancer-${label}`}
                type="button"
                onClick={() => onPresetSelect({ category: selectedCancerType, label })}
                className={`${presetClassName} bg-emerald-50 text-emerald-800 ring-emerald-200 hover:bg-emerald-100`}
              >
                {label}
              </button>
            ))}
          </>
        ) : (
          <span className="text-slate-500">No cancer-specific presets for the current Dx</span>
        )}
        {globalGeneListLabels.length > 0 && <span className="mx-1 h-5 w-px bg-slate-200" aria-hidden="true" />}
        {globalGeneListLabels.map((label) => (
          <button
            key={`global-${label}`}
            type="button"
            onClick={() => onPresetSelect({ category: GLOBAL_GENE_LIST_CATEGORY, label })}
            className={`${presetClassName} bg-slate-50 text-slate-700 ring-slate-200 hover:bg-slate-100`}
          >
            {label}
          </button>
        ))}
        <span className="flex-1" />
        <label className="inline-flex items-center gap-2 text-slate-600">
          <input
            type="checkbox"
            checked={enableCatVrsQueries}
            onChange={(event) => onEnableCatVrsQueriesChange(event.target.checked)}
          />
          Cat-VRS queries
        </label>
      </div>
    </div>
  );
}
