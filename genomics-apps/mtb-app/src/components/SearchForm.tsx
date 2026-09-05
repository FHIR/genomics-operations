import { useState } from 'react';

interface SearchFormProps {
  searchInput: string;
  setSearchInput: (value: string) => void;
  handleSearch: () => void;
  subjectId: string;
  setSubjectId: (value: string) => void;
}

const patientIDs: string[] = [
    "ABC123","NA19240","ABC456",
    "NA18498", "NA19247", "m123",
    "ABC789", "NA18499", "NA19256",
    "CA12345", "NA18870", "NB6TK328",
    "HCC1143", "NA18871", "NB6TK329",
    "HG00403", "NA19190", "XYZ123",
    "HG00406", "NA19210", "XYZ234",
    "HG02657", "NA19238", "XYZ345",
    "huC30902", "NA19239", "L2345"
].sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));

export default function SearchForm({ searchInput, setSearchInput, handleSearch, subjectId, setSubjectId }: SearchFormProps) {
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
    <div className="bg-gray-100 p-8 rounded-lg mb-8 shadow-sm">
      {/* MRN Dropdown field */}
      <div className="mb-4">
        <label htmlFor="mrn-select" className="block text-lg mb-2">MRN (Medical Record Number)</label>
        {!isCustomMRN ? (
          <div className="flex gap-2">
            <select
              id="mrn-select"
              value={subjectId}
              onChange={(e) => handleMRNChange(e.target.value)}
              className="flex-1 bg-white text-black p-4 rounded-lg text-lg border border-gray-200"
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
              className="flex-1 bg-white text-black p-4 rounded-lg text-lg border border-gray-200"
            />
            <button
              onClick={() => {
                setIsCustomMRN(false);
                setSubjectId('');
                setCustomMRN('');
              }}
              className="bg-gray-500 hover:bg-gray-600 text-white px-4 py-2 rounded-lg text-sm"
            >
              Back to List
            </button>
          </div>
        )}
      </div>

      {/* Gene Symbols section */}
      <p className="text-sm text-gray-600 mb-4">Enter Gene Symbols or Genomic Ranges (comma-separated). Ranges must be in zero-based RefSeq:Integer-range format (e.g. &apos;NC_000007.14:55019016-55211628&apos;)</p>
      <div className="flex gap-4">
        {/* Search input field */}
        <input
          type="text"
          value={searchInput}
          onChange={(e) => setSearchInput(e.target.value)}
          className="flex-1 bg-white text-black p-4 rounded-lg text-lg border border-gray-200"
          onKeyDown={(e) => e.key === 'Enter' && handleSearch()} // Trigger search on Enter key
          placeholder="e.g., BRAF, EGFR, ALK"
        />
        {/* Search button */}
        <button
          onClick={handleSearch}
          className="bg-blue-600 hover:bg-blue-700 text-white px-8 py-4 rounded-lg text-lg"
        >
          Search
        </button>
      </div>
      <div className="mt-2 text-sm text-gray-600">
        <p>Examples:</p>
        <ul className="list-disc pl-5 space-y-1">
          <li><span className="font-mono bg-gray-200 px-1 rounded">NC_000007.14:140713327-140924929</span> - BRAF</li>
          <li><span className="font-mono bg-gray-200 px-1 rounded">NC_000007.14:55174721-55174820</span> - EGFR exon 19</li>
          <li><span className="font-mono bg-gray-200 px-1 rounded">NC_000002.12:29190991-29921589</span> - ALK</li>
          <li><span className="font-mono bg-gray-200 px-1 rounded">BRAF, EGFR, ALK</span> - Search multiple genes at once</li>
        </ul>
      </div>
    </div>
  );
}