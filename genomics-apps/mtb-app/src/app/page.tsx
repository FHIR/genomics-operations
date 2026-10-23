'use client';

import Link from "next/link"
import { useEffect, useLayoutEffect, useState, useRef, useCallback, useMemo } from 'react';
import Image from 'next/image';
import { Variant } from '@/types/variants';
import { findSubjectVariantsWithCache } from '@/services/cachedVariantService';
import SearchForm from '../components/SearchForm';
import PatientBanner from '../components/PatientBanner';
import CaseSummary from '../components/caseSummary/CaseSummary';
import SearchStatus from '../components/SearchStatus';
import ResultsTable, { HighlightRequest } from '../components/ResultsTable';
import SelectionTray from '../components/SelectionTray';
import PathwaysView from '../components/pathways/PathwaysView';
import { loadPathways } from '@/services/pathwayService';
import { CooccurrenceProfile, findCooccurrenceMatches, groupMatchesByVariant, loadCooccurrenceProfiles, restrictMatchesToVariants } from '@/services/cooccurrenceService';
import { getVariantGene, getVariantShortLabel } from '@/lib/variantDisplay';
import { CaseDx, CrosswalkRow, getGenericDxOptions, loadTumorTypeCrosswalk } from "@/cancerFilter/tumorTypeCrosswalk";
import { useClinicalRecord } from "@/lib/useClinicalRecord";
import { getPatientDxOptions, getPrimaryCancerConditions } from "@/lib/clinicalSummary";
import RegionLoader from "@/cancerFilter/regionLoader";
import { PresetSelection } from "@/cancerFilter/geneListKb";
import FilterSidebar from "@/cancerFilter/FilterSidebar";
import { DEFAULT_ACTIONABILITY_FILTER, FilterCriteria } from "@/cancerFilter/FilterSidebar";
import { applyFiltersToVariants } from "@/cancerFilter/filterUtils";
import { isPhenotypeMatchForImplication, isALevelEvidence } from "@/cancerFilter/phenotypeUtils";
import { HOW_TO_USE_GROUPS, HOW_TO_USE_TITLE, HowToUseLinkItem } from './howToUseContent';
// import EmailSubscription from "@/components/EmailSubscription";
// import FeedbackForm from "@/components/FeedbackForm";

function isHowToUseLinkItem(item: string | HowToUseLinkItem): item is HowToUseLinkItem {
  return typeof item !== 'string';
}

type CaseView = 'summary' | 'genomic' | 'synthesis';

const CASE_VIEWS: { view: CaseView; label: string; disabledReason?: string }[] = [
  { view: 'summary', label: 'Case Summary' },
  { view: 'genomic', label: 'Genomic Findings' },
  { view: 'synthesis', label: 'Synthesis', disabledReason: 'Not built yet' },
];

function hasActiveSidebarFilters(filters: FilterCriteria, cancerType: string) {
  const actionabilityIsActive = filters.actionableFilter !== 'None' &&
    (filters.actionableFilter !== DEFAULT_ACTIONABILITY_FILTER || Boolean(cancerType));

  return actionabilityIsActive ||
    filters.selectedImpacts.length > 0 ||
    filters.selectedConsequences.length > 0 ||
    filters.selectedDxSignificances.length > 0 ||
    filters.selectedDxStars.length > 0 ||
    filters.selectedTxEvidence.length > 0 ||
    filters.selectedTxMedications.length > 0 ||
    filters.selectedTxImplications.length > 0 ||
    filters.selectedTxPhenotypes.length > 0;
}

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
  const [enableCatVrsQueries, setEnableCatVrsQueries] = useState(false);
  const [isHowToUseOpen, setIsHowToUseOpen] = useState(false);
  // The case Dx is the source of truth for tumor type: it sets the gene-list presets, implication
  // filtering and the oncogenicity predictor's tumorType
  const [caseDx, setCaseDx] = useState<CaseDx | null>(null);
  const [tumorTypeCrosswalk, setTumorTypeCrosswalk] = useState<CrosswalkRow[]>([]);
  const [selectedPreset, setSelectedPreset] = useState<PresetSelection | null>(null);
  const [presetRequestId, setPresetRequestId] = useState(0);
  const [isSearchWorkspaceCollapsed, setIsSearchWorkspaceCollapsed] = useState(false);
  const [isFilterOpen, setIsFilterOpen] = useState(false);
  const [activeView, setActiveView] = useState<'results' | 'pathways'>('results');
  // The case views below the banner; Genomic findings is the default, so the app opens as it always has
  const [caseView, setCaseView] = useState<CaseView>('genomic');
  const [selectedVariantIds, setSelectedVariantIds] = useState<Set<string>>(new Set());
  const [visibleVariantIds, setVisibleVariantIds] = useState<string[]>([]);
  const [highlightRequest, setHighlightRequest] = useState<HighlightRequest | null>(null);
  const [trayMessage, setTrayMessage] = useState<string | null>(null);
  const [cooccurrenceProfiles, setCooccurrenceProfiles] = useState<CooccurrenceProfile[]>([]);
  const [currentFilters, setCurrentFilters] = useState<FilterCriteria>({
    selectedImpacts: [],
    selectedConsequences: [],
    selectedDxSignificances: [],
    selectedDxStars: [],
    selectedTxEvidence: [],
    selectedTxMedications: [],
    selectedTxImplications: [],
    selectedTxPhenotypes: [],
    actionableFilter: DEFAULT_ACTIONABILITY_FILTER,
  });

  const selectedCancerType = caseDx?.geneListCategory ?? '';
  const { state: clinicalState, retry: retryClinicalRecord } = useClinicalRecord(subjectId);
  const clinicalRecord = clinicalState.status === 'loaded' ? clinicalState.record : null;

  useEffect(() => {
    loadTumorTypeCrosswalk()
      .then(setTumorTypeCrosswalk)
      .catch((error) => console.error('Error loading tumor type crosswalk:', error));
  }, []);

  const primaryCancerConditions = useMemo(
    () => (clinicalRecord ? getPrimaryCancerConditions(clinicalRecord.resources, tumorTypeCrosswalk) : []),
    [clinicalRecord, tumorTypeCrosswalk]
  );
  const patientDxOptions = useMemo(
    () => getPatientDxOptions(primaryCancerConditions, tumorTypeCrosswalk),
    [primaryCancerConditions, tumorTypeCrosswalk]
  );
  const genericDxOptions = useMemo(() => getGenericDxOptions(tumorTypeCrosswalk), [tumorTypeCrosswalk]);
  // Stage and disease status follow the selected patient cancer, or the patient's main one after an override
  const bannerCondition = primaryCancerConditions.find((condition) => condition.id === caseDx?.conditionId)
    ?? primaryCancerConditions[0];

  // A new patient's Dx starts as their own primary cancer, or empty when they have no clinical data
  useEffect(() => {
    setCaseDx(patientDxOptions[0] ?? null);
  }, [patientDxOptions]);

  // Load the pathway diagrams in the background so the Pathways view opens without a loading step
  useEffect(() => {
    loadPathways().catch((error) => console.error('Error loading pathway diagrams:', error));
  }, []);

  useEffect(() => {
    loadCooccurrenceProfiles()
      .then(setCooccurrenceProfiles)
      .catch((error) => console.error('Error loading co-occurrence knowledge base:', error));
  }, []);

  // Co-occurrence is checked against all results, so filters never hide it; it reruns as results arrive
  const cooccurrenceMatches = useMemo(
    () => findCooccurrenceMatches(results, cooccurrenceProfiles),
    [results, cooccurrenceProfiles]
  );
  const cooccurrencesByVariant = useMemo(() => groupMatchesByVariant(cooccurrenceMatches), [cooccurrenceMatches]);

  const selectedVariants = useMemo(
    () => results.filter((variant) => variant.id && selectedVariantIds.has(variant.id)),
    [results, selectedVariantIds]
  );
  const visibleVariantIdSet = useMemo(() => new Set(visibleVariantIds), [visibleVariantIds]);
  const hiddenSelectedCount = selectedVariants.filter((variant) => !visibleVariantIdSet.has(variant.id as string)).length;
  const selectedCooccurrenceMatches = useMemo(
    () => restrictMatchesToVariants(cooccurrenceMatches, selectedVariantIds),
    [cooccurrenceMatches, selectedVariantIds]
  );

  // Keep the tab bar where it is on screen when switching views (or at the top of the window if it
  // was scrolled out of view), and return to the previous scroll position in the results table
  const tabBarRef = useRef<HTMLDivElement | null>(null);
  const pendingTabBarTopRef = useRef<number | null>(null);
  const resultsScrollYRef = useRef<number | null>(null);

  const switchView = useCallback((view: 'results' | 'pathways') => {
    if (view === activeView) {
      return;
    }
    if (activeView === 'results') {
      resultsScrollYRef.current = window.scrollY;
    }
    pendingTabBarTopRef.current = tabBarRef.current?.getBoundingClientRect().top ?? null;
    setActiveView(view);
  }, [activeView]);

  useLayoutEffect(() => {
    const previousTop = pendingTabBarTopRef.current;
    pendingTabBarTopRef.current = null;
    if (previousTop === null || !tabBarRef.current) {
      return;
    }
    if (activeView === 'results' && resultsScrollYRef.current !== null) {
      window.scrollTo({ top: resultsScrollYRef.current });
      return;
    }
    window.scrollBy({ top: tabBarRef.current.getBoundingClientRect().top - Math.max(previousTop, 0) });
  }, [activeView]);

  const toggleVariantSelection = useCallback((variantId: string) => {
    setSelectedVariantIds((current) => {
      const next = new Set(current);
      if (next.has(variantId)) {
        next.delete(variantId);
      } else {
        next.add(variantId);
      }
      return next;
    });
  }, []);

  const selectVariants = useCallback((variantIds: string[]) => {
    setSelectedVariantIds((current) => new Set([...current, ...variantIds]));
  }, []);

  const clearSelection = useCallback(() => {
    setSelectedVariantIds(new Set());
    setActiveView('results');
  }, []);

  const showInResults = useCallback((variantId: string) => {
    setTrayMessage(null);
    setActiveView('results');
    setHighlightRequest((current) => ({ variantId, requestId: (current?.requestId ?? 0) + 1 }));
  }, []);

  const handleHighlightMissing = useCallback((variantId: string) => {
    const variant = results.find((candidate) => candidate.id === variantId);
    const label = variant ? `${getVariantGene(variant) ?? variant.range} ${getVariantShortLabel(variant)}` : 'That variant';
    setTrayMessage(`${label} is hidden by current filters.`);
  }, [results]);

  useEffect(() => {
    if (!trayMessage) {
      return;
    }
    const timer = window.setTimeout(() => setTrayMessage(null), 5000);
    return () => window.clearTimeout(timer);
  }, [trayMessage]);

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
    if (!isHowToUseOpen) {
      return;
    }

    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        setIsHowToUseOpen(false);
      }
    };

    document.addEventListener('keydown', handleKeyDown);

    return () => {
      document.removeEventListener('keydown', handleKeyDown);
    };
  }, [isHowToUseOpen]);

  useEffect(() => {
    setSelectedPreset(null);
    setPresetRequestId(0);
  }, [selectedCancerType]);

  // A different patient starts with no genomic results; the search terms are kept so they can be rerun
  const handleSubjectIdChange = useCallback((value: string) => {
    abortControllersRef.current.forEach((controller) => controller.abort());
    abortControllersRef.current.clear();
    setResults([]);
    setFilteredResults([]);
    setSearchStatus({});
    setInvalidRanges([]);
    setSelectedVariantIds(new Set());
    setActiveView('results');
    setSubjectId(value);
  }, []);

  const handlePresetSelect = useCallback((selection: PresetSelection) => {
    setSelectedPreset(selection);
    setPresetRequestId((currentRequestId) => currentRequestId + 1);
  }, []);

  const handleSearchInputChange = useCallback((value: string) => {
    setSelectedPreset(null);
    setPresetRequestId(0);
    setSearchInput(value);
  }, []);

  // Unified filter function that combines actionability and sidebar filters.
  const applyUnifiedFilters = useCallback((
    variants: Variant[],
    filters: FilterCriteria,
    cancerType: string,
  ): Variant[] => {
    try {
      let filtered = applyFiltersToVariants(variants, filters);
      const actionableLabel = filters.actionableFilter;

      if (actionableLabel === "Actionable, this tumor type") {
        if (!cancerType) {
          return filtered;
        }

        filtered = filtered.map(variant => ({
          ...variant,
          txImplications: variant.txImplications?.filter(tx => {
            const hasALevel = isALevelEvidence(tx.evidenceLevel);
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
      } else if (actionableLabel === "Actionable, any tumor type") {
        filtered = filtered.map(variant => ({
          ...variant,
          txImplications: variant.txImplications?.filter(tx => {
            const hasALevel = isALevelEvidence(tx.evidenceLevel);

            return hasALevel;
          }),
          dxImplications: variant.dxImplications?.filter(dx => {
            const hasALevel = isALevelEvidence(dx.evidenceLevel);

            return hasALevel;
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
    const filtered = applyUnifiedFilters(results, filters, selectedCancerType);
    setFilteredResults(filtered);
  }, [results, selectedCancerType, applyUnifiedFilters]);

  // Apply filters whenever currentFilters, results, or related dependencies change
  useEffect(() => {
    if (results.length > 0) {
      const filtered = applyUnifiedFilters(results, currentFilters, selectedCancerType);
      setFilteredResults(filtered);
    } else {
      setFilteredResults([]);
    }
  }, [results, currentFilters, selectedCancerType, applyUnifiedFilters]);

  const hasActiveFilters = hasActiveSidebarFilters(currentFilters, selectedCancerType);
  const trimmedSearchInput = searchInput.trim();
  const searchTermCount = trimmedSearchInput
    ? trimmedSearchInput.split(',').map(term => term.trim()).filter(Boolean).length
    : 0;
  const searchPreview = trimmedSearchInput
    ? (trimmedSearchInput.length > 60 ? `${trimmedSearchInput.slice(0, 60)}...` : trimmedSearchInput)
    : 'No search terms loaded';
  const hasSearchActivity = Object.keys(searchStatus).length > 0 || invalidRanges.length > 0;

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
    setSelectedVariantIds(new Set());
    setActiveView('results');

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
        const processedVariants = await findSubjectVariantsWithCache(term, subjectId, searchId, enableCatVrsQueries);

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
    <div className="min-h-screen bg-white px-4 py-6 text-black sm:px-6 lg:px-8">
      {/* One centered container that fills the window up to 1400px; sections adapt to its width */}
      <main className={`mx-auto w-full max-w-[1400px] ${selectedVariants.length > 0 ? 'pb-24' : ''}`}>
        <div className="mb-8 flex flex-wrap items-center justify-between gap-4">
          <div className="flex flex-wrap items-center gap-4">
            <h1 className="text-4xl font-bold sm:text-5xl">Molecular Tumor Board</h1>
            <button
              type="button"
              onClick={() => setIsHowToUseOpen(true)}
              className="inline-flex items-center rounded-md border border-blue-200 bg-blue-50 px-3 py-1.5 text-sm font-semibold text-blue-700 transition-colors hover:bg-blue-100 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
            >
              How to use
            </button>
          </div>
          <Link href="https://elimu.io/" target="_blank" rel="noopener noreferrer">
            <Image
              src="/elimu.png"
              alt="Elimu"
              width={200}
              height={200}
              // Smaller on narrow windows so it stays beside the title
              className="h-auto w-32 cursor-pointer lg:w-[200px]"
            />
          </Link>
        </div>

        <PatientBanner
          subjectId={subjectId}
          onSubjectIdChange={handleSubjectIdChange}
          patientOptions={patientDxOptions}
          genericOptions={genericDxOptions}
          selectedDx={caseDx}
          onDxChange={setCaseDx}
          clinicalState={clinicalState}
          onRetry={retryClinicalRecord}
          bannerCondition={bannerCondition}
        />

        <div className="mb-6 flex flex-wrap gap-2" role="tablist" aria-label="Case views">
          {CASE_VIEWS.map(({ view, label, disabledReason }) => {
            const isActive = caseView === view;
            return (
              <button
                key={view}
                type="button"
                role="tab"
                aria-selected={isActive}
                disabled={Boolean(disabledReason)}
                title={disabledReason}
                onClick={() => setCaseView(view)}
                className={`rounded-lg border px-5 py-2 text-sm font-semibold ${isActive ? 'border-slate-900 bg-slate-900 text-white' : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'} disabled:cursor-not-allowed disabled:opacity-50`}
              >
                {label}
              </button>
            );
          })}
        </div>

        {caseView === 'summary' && (
          // White cards on a light gray ground give the summary some depth
          <div className="mb-8 rounded-3xl bg-slate-100 p-3 ring-1 ring-inset ring-slate-200 sm:p-5">
            <CaseSummary
              subjectId={subjectId}
              clinicalState={clinicalState}
              onRetry={retryClinicalRecord}
              crosswalk={tumorTypeCrosswalk}
            />
          </div>
        )}

        <RegionLoader
          presetSelection={selectedPreset}
          requestId={presetRequestId}
          onRegionsLoaded={(regions) => {
            setSearchInput(regions.join(", "));
            setSelectedPreset(null);
            setPresetRequestId(0);
          }}
        />

        {/* Genomic findings stay mounted on other tabs, so searches, results and predictions survive switching */}
        <div className={caseView === 'genomic' ? '' : 'hidden'}>
          <div className="mb-8 overflow-hidden rounded-2xl border-2 border-slate-400 bg-gradient-to-br from-slate-100 to-slate-50 shadow-sm">
            <div className="flex flex-wrap items-center justify-between gap-3 border-b-2 border-slate-300 px-5 py-4">
              <div>
                <h2 className="text-lg font-semibold text-gray-900">Search Workspace</h2>
                <p className="mt-1 text-sm text-gray-600">
                  {selectedCancerType ? `Cancer type: ${selectedCancerType}` : 'Select a Dx with gene lists to load cancer-specific preset search terms.'}
                  {' '}
                  {searchTermCount > 0 ? `• ${searchTermCount} search term${searchTermCount === 1 ? '' : 's'}` : ''}
                </p>
              </div>
              <button
                type="button"
                onClick={() => setIsSearchWorkspaceCollapsed((currentValue) => !currentValue)}
                className="inline-flex items-center gap-2 rounded-md border border-blue-700 bg-blue-600 px-4 py-2.5 text-sm font-semibold text-white shadow-[0_3px_0_0_rgb(29_78_216)] transition-[transform,box-shadow,background-color] hover:bg-blue-700 hover:shadow-[0_2px_0_0_rgb(30_64_175)] active:translate-y-px active:shadow-[0_1px_0_0_rgb(30_64_175)] focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
              >
                {isSearchWorkspaceCollapsed ? 'Show Search Panel' : 'Hide Search Panel'}
              </button>
            </div>

            {isSearchWorkspaceCollapsed ? (
              <div className="flex flex-wrap items-center justify-between gap-3 px-5 py-4 text-sm text-gray-600">
                <span className="min-w-0 flex-1 truncate">{searchPreview}</span>
                <div className="flex items-center gap-3 text-xs font-medium uppercase tracking-wide text-gray-500">
                  <span>{hasSearchActivity ? 'Search status available' : 'No search activity yet'}</span>
                  <span>{hasActiveFilters ? 'Filters active' : 'No active filters'}</span>
                </div>
              </div>
            ) : (
              <div className="grid grid-cols-1 items-start gap-0 xl:grid-cols-[minmax(0,1.65fr)_14rem]">
                <SearchForm
                  searchInput={searchInput}
                  setSearchInput={handleSearchInputChange}
                  handleSearch={handleSearch}
                  selectedCancerType={selectedCancerType}
                  onPresetSelect={handlePresetSelect}
                  enableCatVrsQueries={enableCatVrsQueries}
                  onEnableCatVrsQueriesChange={setEnableCatVrsQueries}
                  embedded={true}
                  className="p-5"
                />
                <div className="border-t-2 border-slate-300 xl:border-l-2 xl:border-t-0">
                  <SearchStatus
                    searchStatus={searchStatus}
                    invalidRanges={invalidRanges}
                    onCancelSearch={handleCancelSearch}
                    embedded={true}
                    className="max-h-[32rem]"
                  />
                </div>
              </div>
            )}
          </div>
          <FilterSidebar
            onFilterChange={handleFilterChange}
            isOpen={isFilterOpen}
            onClose={() => setIsFilterOpen(false)}
            results={results}
            currentFilters={currentFilters}
          />
          <div ref={tabBarRef} className="flex flex-wrap gap-1 border-b border-slate-300" role="tablist">
            {(['results', 'pathways'] as const).map((view) => {
              const isActive = activeView === view;
              const isDisabled = view === 'pathways' && selectedVariants.length === 0;
              return (
                <button
                  key={view}
                  type="button"
                  role="tab"
                  aria-selected={isActive}
                  disabled={isDisabled}
                  title={isDisabled ? 'Select variants in the results table first' : undefined}
                  onClick={() => switchView(view)}
                  className={`-mb-px rounded-t-lg border px-5 py-2.5 text-sm font-semibold ${isActive ? 'border-slate-300 border-b-gray-100 bg-gray-100 text-gray-900' : 'border-transparent text-gray-500 hover:text-gray-800'} disabled:cursor-not-allowed disabled:opacity-50`}
                >
                  {view === 'results' ? `Results (${filteredResults.length})` : `Pathways${selectedVariants.length ? ` (${selectedVariants.length} selected)` : ''}`}
                </button>
              );
            })}
          </div>
          {/* The results table stays mounted so its state (e.g., oncogenicity results) survives switching tabs */}
          <div className={activeView === 'results' ? '' : 'hidden'}>
            <ResultsTable
              results={filteredResults}
              predictorTumorType={caseDx?.predictorTumorType}
              onToggleFilters={() => setIsFilterOpen(!isFilterOpen)}
              hasActiveFilters={hasActiveFilters}
              enableCatVrsQueries={enableCatVrsQueries}
              onEnableCatVrsQueriesChange={setEnableCatVrsQueries}
              cooccurrencesByVariant={cooccurrencesByVariant}
              selectedVariantIds={selectedVariantIds}
              onToggleVariantSelection={toggleVariantSelection}
              onSelectVariants={selectVariants}
              onVisibleVariantIdsChange={setVisibleVariantIds}
              highlightRequest={highlightRequest}
              onHighlightMissing={handleHighlightMissing}
            />
          </div>
          {activeView === 'pathways' && (
            // At least a window tall, so switching from the long results table doesn't pull the page up
            <div className="min-h-screen rounded-2xl border-2 border-slate-400 bg-gray-100 p-3 shadow-sm sm:p-6">
              <PathwaysView
                selectedVariants={selectedVariants}
                cooccurrenceMatches={selectedCooccurrenceMatches}
                cooccurrencesByVariant={cooccurrencesByVariant}
                onShowInResults={showInResults}
              />
            </div>
          )}
        </div>
        {caseView === 'genomic' && activeView === 'results' && (
          <SelectionTray
            selectedVariants={selectedVariants}
            hiddenCount={hiddenSelectedCount}
            message={trayMessage}
            onClear={clearSelection}
            onViewPathways={() => switchView('pathways')}
          />
        )}
        {isHowToUseOpen && (
          <div
            className="fixed inset-0 z-50 flex items-center justify-center bg-black/45 px-4 py-8"
            onClick={() => setIsHowToUseOpen(false)}
          >
            <div
              role="dialog"
              aria-modal="true"
              aria-labelledby="how-to-use-title"
              className="max-h-[85vh] w-full max-w-3xl overflow-y-auto rounded-2xl bg-white p-6 shadow-2xl"
              onClick={(event) => event.stopPropagation()}
            >
              <div className="mb-6 flex items-start justify-between gap-4">
                <div>
                  <h2 id="how-to-use-title" className="text-2xl font-bold text-gray-900">
                    {HOW_TO_USE_TITLE}
                  </h2>
                  <p className="mt-2 text-sm text-gray-600">
                    These notes summarize the current search, filtering, and results behavior in the app.
                  </p>
                </div>
                <button
                  type="button"
                  onClick={() => setIsHowToUseOpen(false)}
                  className="rounded-md border border-gray-300 px-3 py-1.5 text-sm font-medium text-gray-700 transition-colors hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
                >
                  Close
                </button>
              </div>
              <div className="space-y-8">
                {HOW_TO_USE_GROUPS.map((group, groupIndex) => (
                  // Each tab's guidance is a group under its own heading; the untitled first group applies to all tabs
                  <div key={group.title ?? `group-${groupIndex}`} className="space-y-6">
                    {group.title && (
                      <h3 className="border-b-2 border-slate-200 pb-1.5 text-xl font-bold text-gray-900">{group.title}</h3>
                    )}
                    {group.sections.map((section) => (
                      <section key={`${group.title ?? 'general'}-${section.title}`}>
                        <h4 className="text-lg font-semibold text-gray-900">{section.title}</h4>
                        <ul className="mt-3 list-disc space-y-2 pl-5 text-sm leading-6 text-gray-700">
                          {section.items.map((item) => (
                            <li key={typeof item === 'string' ? item : `${section.title}-${item.linkLabel}`}>
                              {isHowToUseLinkItem(item) ? (
                                <>
                                  {item.textBefore}
                                  <a
                                    href={item.linkUrl}
                                    target="_blank"
                                    rel="noopener noreferrer"
                                    className="font-medium text-blue-600 hover:text-blue-800 hover:underline"
                                  >
                                    {item.linkLabel}
                                  </a>
                                  {item.textAfter}
                                </>
                              ) : item}
                            </li>
                          ))}
                        </ul>
                      </section>
                    ))}
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}
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
