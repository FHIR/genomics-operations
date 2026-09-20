interface SearchStatusProps {
  searchStatus: Record<string, string>;
  invalidRanges: string[];
  onCancelSearch: () => void;
  embedded?: boolean;
  className?: string;
}

export default function SearchStatus({
  searchStatus,
  invalidRanges,
  onCancelSearch,
  embedded = false,
  className = '',
}: SearchStatusProps) {
  // Determine if any range is currently searching
  const isSearching = Object.values(searchStatus).includes('searching');
  const hasSearchStatus = Object.keys(searchStatus).length > 0;
  const hasInvalidRanges = invalidRanges.length > 0;

  if (!hasSearchStatus && !hasInvalidRanges && !embedded) {
    return null;
  }

  return (
    <>
      {/* Search Status Section */}
      {hasSearchStatus ? (
        <div className={`${embedded ? '' : 'bg-gray-100 rounded-lg shadow-sm'} p-4 overflow-hidden ${className}`.trim()}>
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
      ) : embedded ? (
        <div className={`p-4 ${className}`.trim()}>
          <h3 className="mb-2 text-base font-medium">Search Status</h3>
          <p className="text-sm text-gray-500">Status updates will appear here after you run a search.</p>
        </div>
      ) : null}

      {/* Invalid Ranges Section */}
      {hasInvalidRanges && (
        <div className={`${embedded ? 'mt-4 border border-red-200 bg-red-50' : 'mb-4 border border-red-200 bg-red-50 shadow-sm'} rounded-lg p-4`}>
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
