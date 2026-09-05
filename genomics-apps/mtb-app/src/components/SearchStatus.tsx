interface SearchStatusProps {
  searchStatus: Record<string, string>;
  invalidRanges: string[];
  onCancelSearch: () => void;
  className?: string;
}

export default function SearchStatus({
  searchStatus,
  invalidRanges,
  onCancelSearch,
  className = '',
}: SearchStatusProps) {
  // Determine if any range is currently searching
  const isSearching = Object.values(searchStatus).includes('searching');

  return (
    <>
      {/* Search Status Section */}
      {Object.keys(searchStatus).length > 0 && (
        <div className={`bg-gray-100 p-4 rounded-lg shadow-sm overflow-hidden ${className}`.trim()}>
          <h3 className="mb-2 text-base font-medium">Search Status</h3>
          <div className="space-y-1.5 overflow-y-auto pr-1 max-h-[24rem] text-sm">
            {Object.entries(searchStatus).map(([range, status]) => (
              <div key={range} className="flex items-start leading-5">
                <span className="mr-2 font-medium break-all">{range}:</span>
                {status === 'searching' && (
                  <span className="text-yellow-600">Searching...</span>
                )}
                {status === 'found' && (
                  <span className="text-green-600">Results found</span>
                )}
                {status === 'no results' && (
                  <span className="text-gray-600">No results found</span>
                )}
                {status === 'error' && (
                  <span className="text-red-600">Error processing</span>
                )}
                {status === 'invalid_mrn' && (
                  <span className="text-red-600">Invalid MRN - Patient not found</span>
                )}
                {status === 'cancelled' && (
                  <span className="text-orange-600">Cancelled</span>
                )}
              </div>
            ))}
          </div>

          {isSearching && (
            <div className="mt-4">
              <button
                className="rounded bg-red-500 px-3 py-1.5 text-sm text-white transition-colors hover:bg-red-600"
                onClick={onCancelSearch}
              >
                Cancel Search
              </button>
            </div>
          )}
        </div>
      )}

      {/* Invalid Ranges Section */}
      {invalidRanges.length > 0 && (
        <div className="bg-red-50 p-4 rounded-lg mb-4 shadow-sm border border-red-200">
          <h3 className="text-lg font-medium text-red-700 mb-2">Unprocessable Ranges</h3>
          <ul className="list-disc pl-5 text-red-600">
            {invalidRanges.map((range, index) => (
              <li key={index}>{range} - Unrecognized gene or malformed chromosomal range</li>
            ))}
          </ul>
        </div>
      )}
    </>
  );
}
