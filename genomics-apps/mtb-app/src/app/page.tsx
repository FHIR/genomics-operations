'use client';

import Link from "next/link"
import { useEffect, useState, useRef, useCallback } from 'react';
import Image from 'next/image';
import { Variant } from '@/types/variants';
import { findSubjectVariantsWithCache } from '@/services/cachedVariantService';
import SearchForm from '../components/SearchForm';
import SearchStatus from '../components/SearchStatus';
import ResultsTable from '../components/ResultsTable';
import CancerSelect from "@/cancerFilter/cancerSelect";
import ActionableCheckBoxes from "@/cancerFilter/actionableCheckBoxes";
import RegionLoader from "@/cancerFilter/regionLoader";
import FilterSidebar from "@/cancerFilter/FilterSidebar";
import { FilterCriteria } from "@/cancerFilter/FilterSidebar";
import { applyFiltersToVariants } from "@/cancerFilter/filterUtils";
import { isPhenotypeMatchForImplication, isOtherTumorType, isALevelEvidence } from "@/cancerFilter/phenotypeUtils";
// import EmailSubscription from "@/components/EmailSubscription";
// import FeedbackForm from "@/components/FeedbackForm";
export default function Home() {
  // useRef for better search ID generation
  const searchIdRef = useRef(0);
  // useRef for abort controllers to avoid dependency issues
  const abortControllersRef = useRef<Map<string, AbortController>>(new Map());

  // State management using React hooks
  const [searchInput, setSearchInput] = useState('');
  const [subjectId, setSubjectId] = useState('L2345');
  const [results, setResults] = useState<Variant[]>([]);
  const [filteredResults, setFilteredResults] = useState<Variant[]>([]);
  const [searchStatus, setSearchStatus] = useState<Record<string, string>>({});
  const [invalidRanges, setInvalidRanges] = useState<string[]>([]);
  const [selectedCancerType, setSelectedCancerType] = useState("");
  const [selectedLabel, setSelectedLabel] = useState("");
  const [isFilterOpen, setIsFilterOpen] = useState(false);
  const [mtbKbPhenotypes, setMtbKbPhenotypes] = useState<Set<string>>(new Set());
  const [currentFilters, setCurrentFilters] = useState<FilterCriteria>({
    selectedImpacts: [],
    selectedConsequences: [],
    selectedDxSignificances: [],
    selectedDxStars: [],
    selectedTxEvidence: [],
    selectedTxMedications: [],
    selectedTxImplications: [],
    selectedTxPhenotypes: []
  });
  /**
   * Tracks whether regions have been loaded for the current cancer type and label combination.
   * Reset to false when either selection changes to allow reloading.
   */
  const [hasLoadedRegions, setHasLoadedRegions] = useState(false);
  // Cleanup abort controllers when component unmounts
  useEffect(() => {
    const controllers = abortControllersRef.current;
    return () => {
      controllers.forEach(controller => controller.abort());
    };
  }, []);

  const handleCancelSearch = () => {
    // abort all ongoing searches
    abortControllersRef.current.forEach((controller) => {
      controller.abort();
    });
    abortControllersRef.current.clear();

    // update search status for all 'searching' terms to 'cancelled'
    setSearchStatus(prev => {
      const updated = { ...prev };
      Object.keys(updated).forEach(key => {
        if (updated[key] === 'searching') {
          updated[key] = 'cancelled';
        }
      });
      return updated;
    });
  };

  useEffect(() => {
    if (!selectedLabel) {
      setSearchInput("");
    }
  }, [selectedLabel]);

  // Reset regions loading state when either selection changes
  useEffect(() => {
    setHasLoadedRegions(false);
  }, [selectedCancerType, selectedLabel]);

  // Unified filter function that combines both radio and user filters
  const applyUnifiedFilters = useCallback((
    variants: Variant[],
    filters: FilterCriteria,
    cancerType: string,
    actionableLabel: string,
    phenotypes: Set<string>
  ): Variant[] => {
    try {
      // First apply the dynamic filters from FilterSidebar
      let filtered = applyFiltersToVariants(variants, filters);

      // Then apply actionable logic based on radio button selection
      if (actionableLabel === "Actionable, this tumor type") {
        filtered = filtered.map(variant => ({
          ...variant,
          txImplications: variant.txImplications?.filter(tx => {
            // Must have A-level evidence
            const hasALevel = isALevelEvidence(tx.evidenceLevel);
            // Must match cancer type phenotype
            const phenotypeMatch = isPhenotypeMatchForImplication(tx.phenotypicContext, cancerType);

            return hasALevel && phenotypeMatch;
          }),
          dxImplications: variant.dxImplications?.filter(dx => {
            const hasALevel = isALevelEvidence(dx.evidenceLevel);

            const phenotypeMatch = isPhenotypeMatchForImplication(dx.predictedPhenotype, cancerType);

            return hasALevel && phenotypeMatch;
          })
        })).filter(variant =>
          (variant.txImplications && variant.txImplications.length > 0) ||
          (variant.dxImplications && variant.dxImplications.length > 0)
        );
      } else if (actionableLabel === "Actionable, other tumor type") {
        filtered = filtered.map(variant => ({
          ...variant,
          txImplications: variant.txImplications?.filter(tx => {
            const hasALevel = isALevelEvidence(tx.evidenceLevel);

            const isOtherTumor = isOtherTumorType(tx.phenotypicContext, cancerType, phenotypes);

            return hasALevel && isOtherTumor;
          }),
          dxImplications: variant.dxImplications?.filter(dx => {
            const hasALevel = isALevelEvidence(dx.evidenceLevel);

            const isOtherTumor = isOtherTumorType(dx.predictedPhenotype, cancerType, phenotypes);

            return hasALevel && isOtherTumor;
          })
        })).filter(variant =>
          (variant.txImplications && variant.txImplications.length > 0) ||
          (variant.dxImplications && variant.dxImplications.length > 0)
        );
      } else if (actionableLabel === "Possibly actionable") {
        filtered = filtered.map(variant => ({
          ...variant,
          txImplications: variant.txImplications?.filter(tx => {
            // Must NOT have A-level evidence
            const hasALevel = isALevelEvidence(tx.evidenceLevel);

            return !hasALevel;
          }),
          dxImplications: variant.dxImplications?.filter(dx => {
            const hasALevel = isALevelEvidence(dx.evidenceLevel);

            return !hasALevel;
          })
        })).filter(variant =>
          (variant.txImplications && variant.txImplications.length > 0) ||
          (variant.dxImplications && variant.dxImplications.length > 0)
        );
      }

      return filtered;
    } catch (error) {
      console.error('Error applying unified filters:', error);
      return variants; // return original variants on error
    }
  }, []); // Empty dependency array since the function doesn't depend on any external values

  const handleFilterChange = useCallback((filters: FilterCriteria) => {
    setCurrentFilters(filters);
    const filtered = applyUnifiedFilters(results, filters, selectedCancerType, selectedLabel, mtbKbPhenotypes);
    setFilteredResults(filtered);
  }, [results, selectedCancerType, selectedLabel, mtbKbPhenotypes, applyUnifiedFilters]);

  // Update filter state when radio button selection changes
  useEffect(() => {
    setCurrentFilters(prevFilters => {
      const updatedFilters = { ...prevFilters };
      if (selectedLabel === "Actionable, this tumor type" || selectedLabel === "Actionable, other tumor type") {
        // Add A-level evidence filter from radio button
        // We need to match the actual evidence levels from the data
        const radioFilters: { value: string; source: 'radio' }[] = [];

        // Check what A-level evidence levels actually exist in the results
        const existingEvidenceLevels = new Set<string>();
        results.forEach(variant => {
          variant.txImplications?.forEach(tx => {
            if (tx.evidenceLevel && isALevelEvidence(tx.evidenceLevel)) {
              existingEvidenceLevels.add(tx.evidenceLevel);
            }
          });
          variant.dxImplications?.forEach(dx => {
            if (dx.evidenceLevel && isALevelEvidence(dx.evidenceLevel)) {
              existingEvidenceLevels.add(dx.evidenceLevel);
            }
          });
        });

        // Add radio filters for all A-level evidence that actually exists
        existingEvidenceLevels.forEach(evidenceLevel => {
          radioFilters.push({ value: evidenceLevel, source: 'radio' as const });
        });

        // Remove existing radio filters and add the ones that actually exist in data
        updatedFilters.selectedTxEvidence = [
          ...updatedFilters.selectedTxEvidence.filter(f => f.source !== 'radio'),
          ...radioFilters
        ];
        updatedFilters.actionableFilter = selectedLabel;
        updatedFilters.evidenceLevelSource = 'radio';
      } else if (selectedLabel === "Possibly actionable") {
        // For "Possibly actionable", lock all NON-A evidence levels that exist in the data
        const radioFilters: { value: string; source: 'radio' }[] = [];

        // Check what NON-A evidence levels actually exist in the results
        const existingNonAEvidenceLevels = new Set<string>();
        results.forEach(variant => {
          variant.txImplications?.forEach(tx => {
            if (tx.evidenceLevel && !isALevelEvidence(tx.evidenceLevel)) {
              existingNonAEvidenceLevels.add(tx.evidenceLevel);
            }
          });
          variant.dxImplications?.forEach(dx => {
            if (dx.evidenceLevel && !isALevelEvidence(dx.evidenceLevel)) {
              existingNonAEvidenceLevels.add(dx.evidenceLevel);
            }
          });
        });

        // Add radio filters for all non-A evidence levels that actually exist
        existingNonAEvidenceLevels.forEach(evidenceLevel => {
          radioFilters.push({ value: evidenceLevel, source: 'radio' as const });
        });

        // Remove existing radio filters and add the non-A ones
        updatedFilters.selectedTxEvidence = [
          ...updatedFilters.selectedTxEvidence.filter(f => f.source !== 'radio'),
          ...radioFilters
        ];
        updatedFilters.actionableFilter = selectedLabel;
        updatedFilters.evidenceLevelSource = 'radio';
      } else {
        // Clear radio filters when no actionable option is selected
        updatedFilters.selectedTxEvidence = updatedFilters.selectedTxEvidence.filter(f => f.source !== 'radio');
        delete updatedFilters.actionableFilter;
        delete updatedFilters.evidenceLevelSource;
      }
      return updatedFilters;
    });
  }, [selectedLabel, selectedCancerType, results]); // Add results to dependency array

  // Apply filters whenever currentFilters, results, or related dependencies change
  useEffect(() => {
    if (results.length > 0) {
      const filtered = applyUnifiedFilters(results, currentFilters, selectedCancerType, selectedLabel, mtbKbPhenotypes);
      setFilteredResults(filtered);
    } else {
      setFilteredResults([]);
    }
  }, [results, currentFilters, selectedCancerType, selectedLabel, mtbKbPhenotypes, applyUnifiedFilters]);
  const handleSearch = async () => {
    if (!searchInput.trim()) {
      setResults([]);
      setFilteredResults([]);
      setSearchStatus({});
      setInvalidRanges([]);
      return;
    }

    const searchTerms = searchInput
      .split(',')
      .map(term => term.trim())
      .filter(term => term.length > 0)  // Remove empty terms
      .map(term => {
        // Use regex to validate genomic range format: RefSeq:Integer-Integer
        // This pattern matches things like "NC_000007.14:55019016-55211628"
        const rangeRegex = /^[A-Za-z0-9_]+\.[0-9]+:[0-9]+-[0-9]+$/;
        if (rangeRegex.test(term)) {
          return term; // Keep genomic ranges as-is
        }
        return term.toUpperCase(); // Convert gene symbols to uppercase
      });
    // Cancel any ongoing searches
    abortControllersRef.current.forEach(controller => controller.abort());
    abortControllersRef.current.clear();

    // Create new abort controllers for this search
    const newControllers = new Map<string, AbortController>();
    const searchId = ++searchIdRef.current; // Use useRef for better ID generation

    setResults([]);
    setFilteredResults([]);
    setInvalidRanges([]);

    const initialStatus: Record<string, string> = {};
    searchTerms.forEach(term => {
      if (term) {
        initialStatus[term] = 'searching';
        newControllers.set(term, new AbortController());
      }
    });
    setSearchStatus(initialStatus);
    abortControllersRef.current = newControllers;

    // Process each search term sequentially and update results incrementally
    for (const term of searchTerms) {
      if (!term) continue;

      const controller = newControllers.get(term);
      if (!controller) continue;

      try {
        // Use the cached variant service - it handles all implications internally
        const processedVariants = await findSubjectVariantsWithCache(term, subjectId, searchId);

        // Check if this search was cancelled
        if (controller.signal.aborted) {
          continue;
        }

        if (processedVariants.length > 0) {
          // Add all processed variants at once - they already have all implications
          setResults(prev => [...prev, ...processedVariants]);

          if (!controller.signal.aborted) {
            setSearchStatus(prev => ({ ...prev, [term]: 'found' }));
          }
        } else {
          if (!controller.signal.aborted) {
            setSearchStatus(prev => ({ ...prev, [term]: 'no results' }));
          }
        }
      } catch (error) {
        if (!controller.signal.aborted) {
          console.error(`Error searching for ${term}:`, error);

          // Check if the error is related to invalid MRN (400 status)
          const errorMessage = error instanceof Error ? error.message : String(error);
          if (errorMessage.includes('400')) {
            // Check if this might be an MRN issue by seeing if the term looks like a gene
            const isLikelyGene = /^[A-Z0-9]+$/.test(term) || /^NC_\d+\.\d+:\d+-\d+$/.test(term);
            if (isLikelyGene) {
              setSearchStatus(prev => ({ ...prev, [term]: 'invalid_mrn' }));
            } else {
              setSearchStatus(prev => ({ ...prev, [term]: 'error' }));
              setInvalidRanges(prev => [...prev, term]);
            }
          } else {
            setSearchStatus(prev => ({ ...prev, [term]: 'error' }));
            setInvalidRanges(prev => [...prev, term]);
          }
        }
      }
    }
  };
  return (
    <div className="min-h-screen p-8 bg-white text-black">
      <main className="max-w-4xl mx-auto">
        <div className="flex justify-between items-center mb-4">
          <h1 className="text-5xl font-bold">
            Molecular Tumor Board <br /> Genetic Data Viewer
          </h1>
          <Link href="https://elimu.io/" target="_blank" rel="noopener noreferrer">
            <Image
              src="/elimu.png"
              alt="Elimu"
              width={200}
              height={200}
              className="cursor-pointer"
            />
          </Link>
        </div>

        <p className="text-xl text-gray-600 mb-8">Search genetic variants by entering gene symbols or genomic ranges</p>

        <CancerSelect onSelect={setSelectedCancerType} />
        <ActionableCheckBoxes onLabelChange={setSelectedLabel} />

        <RegionLoader
          cancerType={selectedCancerType}
          label={selectedLabel}
          onRegionsLoaded={(regions) => {
            if (!selectedLabel) {
              setSearchInput("");
              return;
            }

            // Only load regions once when both cancer type and label are selected
            // Allow reloading if the user changes the selection
            if (!hasLoadedRegions && regions && regions.length > 0) {
              setSearchInput(regions.join(", "));
              setHasLoadedRegions(true);
            }
          }}
          onPhenotypesLoaded={setMtbKbPhenotypes}
        />

        <div className="mb-8 grid grid-cols-1 items-start gap-6 xl:grid-cols-[minmax(0,1fr)_18rem]">
          <SearchForm
            searchInput={searchInput}
            setSearchInput={setSearchInput}
            handleSearch={handleSearch}
            subjectId={subjectId}
            setSubjectId={setSubjectId}
            className="mb-0"
          />
          <SearchStatus
            searchStatus={searchStatus}
            invalidRanges={invalidRanges}
            onCancelSearch={handleCancelSearch}
            className="max-h-[32rem]"
          />
        </div>
        <FilterSidebar
          onFilterChange={handleFilterChange}
          isOpen={isFilterOpen}
          onClose={() => setIsFilterOpen(false)}
          results={results}
          currentFilters={currentFilters}
        />
        <ResultsTable
          results={filteredResults}
          onToggleFilters={() => setIsFilterOpen(!isFilterOpen)}
        />
        {/*
        <div className="mt-12 grid grid-cols-1 lg:grid-cols-2 gap-8">
          <EmailSubscription />
          <FeedbackForm />
        </div>
        */}
      </main>
    </div>
  );
}
