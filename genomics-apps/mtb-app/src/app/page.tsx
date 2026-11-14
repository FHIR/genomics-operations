'use client';

import Link from "next/link"
import { useEffect, useLayoutEffect, useState, useRef, useCallback, useMemo } from 'react';
import Image from 'next/image';
import { ChevronDown, ChevronUp, ClipboardPlus, Dna, Images, LucideIcon, Network, Sparkles, Table2 } from 'lucide-react';
import { Variant } from '@/types/variants';
import { findSubjectVariantsWithCache } from '@/services/cachedVariantService';
import { PatientNotOnServerError } from '@/services/variantService';
import SearchForm from '../components/SearchForm';
import PatientBanner from '../components/PatientBanner';
import CaseSummary from '../components/caseSummary/CaseSummary';
import SynthesisView from '../components/SynthesisView';
import SearchStatus, { searchSummaryText } from '../components/SearchStatus';
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
import { DEFAULT_ACTIONABILITY_FILTER, effectiveActionabilityFilter, FilterCriteria } from "@/cancerFilter/FilterSidebar";
import { applyFiltersToVariants } from "@/cancerFilter/filterUtils";
import { isALevelEvidence } from "@/cancerFilter/evidenceLevel";
import { loadTumorTypeImplicationCodes, matchesTumorType, TumorTypeCodeIndex } from "@/cancerFilter/tumorTypeCodes";
import { HOW_TO_USE_GROUPS, HOW_TO_USE_TITLE, HowToUseLinkItem } from './howToUseContent';
import { buildReportedVariants, collectReportedVariants } from '@/lib/reportedVariants';
import { layerReportedVariants, needsNormalization, SearchedTerm } from '@/lib/reportedVariantLayering';
import { geneCoordinates, geneForTranscript, NormalizedVariant, normalizeVariant, resolveSearchRange } from '@/services/reportedVariantLookupService';
// import EmailSubscription from "@/components/EmailSubscription";
// import FeedbackForm from "@/components/FeedbackForm";

function isHowToUseLinkItem(item: string | HowToUseLinkItem): item is HowToUseLinkItem {
  return typeof item !== 'string';
}

type CaseView = 'summary' | 'genomic' | 'imaging' | 'synthesis';

// Each tab has an icon in its own color (lighter when the tab is selected); AI Synthesis uses the sparkles icon
// commonly used to mark AI features
// The Genomics Explorer's tooltip says it holds the full data set, unlike the lab's report on the Case Summary
const CASE_VIEWS: { view: CaseView; label: string; icon: LucideIcon; iconClass: string; activeIconClass: string; tooltip?: string }[] = [
  { view: 'summary', label: 'Case Summary', icon: ClipboardPlus, iconClass: 'text-sky-600', activeIconClass: 'text-sky-300' },
  { view: 'genomic', label: 'Genomics Explorer', icon: Dna, iconClass: 'text-emerald-600', activeIconClass: 'text-emerald-300', tooltip: 'Full genomic data set' },
  { view: 'imaging', label: 'Imaging', icon: Images, iconClass: 'text-amber-600', activeIconClass: 'text-amber-300' },
  { view: 'synthesis', label: 'AI Synthesis', icon: Sparkles, iconClass: 'text-violet-600', activeIconClass: 'text-violet-300' },
];

function hasActiveSidebarFilters(filters: FilterCriteria) {
  return filters.actionableFilter !== 'None' ||
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
  // No patient is loaded at launch; the MRN box shows "Select a patient"
  const [subjectId, setSubjectId] = useState('');
  // The genomics server's rows; `results` below adds the patient's reported variants to them
  const [serverResults, setServerResults] = useState<Variant[]>([]);
  const [filteredResults, setFilteredResults] = useState<Variant[]>([]);
  const [searchStatus, setSearchStatus] = useState<Record<string, string>>({});
  const [invalidRanges, setInvalidRanges] = useState<string[]>([]);
  const [enableCatVrsQueries, setEnableCatVrsQueries] = useState(false);
  const [isHowToUseOpen, setIsHowToUseOpen] = useState(false);
  // The case Dx is the source of truth for tumor type: it sets the gene-list presets, implication
  // filtering and the oncogenicity predictor's tumorType
  const [caseDx, setCaseDx] = useState<CaseDx | null>(null);
  const [tumorTypeCrosswalk, setTumorTypeCrosswalk] = useState<CrosswalkRow[]>([]);
  // The cancer-type codes that count for each tumor type, for "Actionable, this tumor type"
  const [tumorTypeCodes, setTumorTypeCodes] = useState<TumorTypeCodeIndex>(new Map());
  const [selectedPreset, setSelectedPreset] = useState<PresetSelection | null>(null);
  const [presetRequestId, setPresetRequestId] = useState(0);
  const [isSearchWorkspaceCollapsed, setIsSearchWorkspaceCollapsed] = useState(false);
  const [isFilterOpen, setIsFilterOpen] = useState(false);
  const [activeView, setActiveView] = useState<'results' | 'pathways'>('results');
  // The case views below the banner. The app opens on Case Summary; once the first patient is chosen it stays there
  // if they have clinical data and moves to Genomics Explorer if they don't (decided once, unless the user picks a tab first); after that the
  // chosen tab stays put when the patient changes
  const [caseView, setCaseView] = useState<CaseView>('summary');
  const initialViewChosenRef = useRef(false);
  const chooseCaseView = useCallback((view: CaseView) => {
    initialViewChosenRef.current = true;
    setCaseView(view);
  }, []);
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

  useEffect(() => {
    if (initialViewChosenRef.current || clinicalState.status === 'idle' || clinicalState.status === 'loading') {
      return;
    }
    initialViewChosenRef.current = true;
    if (clinicalState.status !== 'loaded') {
      setCaseView('genomic');
    }
  }, [clinicalState.status]);
  const clinicalRecord = clinicalState.status === 'loaded' ? clinicalState.record : null;

  useEffect(() => {
    loadTumorTypeCrosswalk()
      .then(setTumorTypeCrosswalk)
      .catch((error) => console.error('Error loading tumor type crosswalk:', error));
    loadTumorTypeImplicationCodes()
      .then(setTumorTypeCodes)
      .catch((error) => console.error('Error loading tumor type implication codes:', error));
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

  // Reported variants: the patient's genetic-report variants join the
  // server's rows for each finished search term whose gene or range covers them, flagged and with credited implications
  const reportedRecords = useMemo(() => {
    if (!clinicalRecord) {
      return [];
    }
    const index = new Map(clinicalRecord.resources.map((resource) => [`${resource.resourceType}/${resource.id}`, resource]));
    return collectReportedVariants(buildReportedVariants(clinicalRecord.resources, index));
  }, [clinicalRecord]);
  const [resolvedRanges, setResolvedRanges] = useState<Map<string, string>>(new Map());
  const [normalizedVariants, setNormalizedVariants] = useState<Map<string, NormalizedVariant | null>>(new Map());
  const [lookedUpGenes, setLookedUpGenes] = useState<Map<string, string>>(new Map());
  const [fusionGeneRanges, setFusionGeneRanges] = useState<Map<string, string[]>>(new Map());
  const requestedLookupsRef = useRef(new Set<string>());
  const searchedTerms: SearchedTerm[] = useMemo(() => Object.entries(searchStatus)
    .filter(([, status]) => status === 'found' || status === 'no results' || status === 'not_on_server')
    .map(([term]) => ({ term, range: resolvedRanges.get(term) })), [searchStatus, resolvedRanges]);

  // Each finished term's range, also for terms the server found nothing in (a reported variant may still be there)
  useEffect(() => {
    if (reportedRecords.length === 0) {
      return;
    }
    searchedTerms.filter(({ term, range }) => !range && !requestedLookupsRef.current.has(`range|${term}`)).forEach(({ term }) => {
      requestedLookupsRef.current.add(`range|${term}`);
      resolveSearchRange(term).then((range) => {
        if (range) {
          setResolvedRanges((previous) => new Map(previous).set(term, range));
        }
      });
    });
  }, [searchedTerms, reportedRecords]);

  // A normalize-variant lookup only for a reported variant the cheap comparisons (SPDI, c.HGVS) can't place
  useEffect(() => {
    for (const record of reportedRecords) {
      const input = record.spdi ?? record.gHgvs ?? record.cHgvs;
      const lookupKey = `normalize|${record.key}`;
      if (!input || requestedLookupsRef.current.has(lookupKey)) {
        continue;
      }
      if (searchedTerms.some((term) => needsNormalization(record, term, serverResults.filter((row) => row.range === term.term)))) {
        requestedLookupsRef.current.add(lookupKey);
        normalizeVariant(input).then((normalized) => setNormalizedVariants((previous) => new Map(previous).set(record.key, normalized)));
      }
    }
  }, [reportedRecords, searchedTerms, serverResults]);

  // A reported fusion's partner coordinates, once a range is searched (a gene search places it by symbol)
  useEffect(() => {
    if (!searchedTerms.some(({ term }) => /^NC_\d+\.\d+:\d+-\d+$/.test(term))) {
      return;
    }
    for (const gene of reportedRecords.flatMap((record) => record.structural?.fusionGenes ?? [])) {
      const lookupKey = `coordinates|${gene}`;
      if (requestedLookupsRef.current.has(lookupKey)) {
        continue;
      }
      requestedLookupsRef.current.add(lookupKey);
      geneCoordinates(gene).then((ranges) => setFusionGeneRanges((previous) => new Map(previous).set(gene, ranges)));
    }
  }, [reportedRecords, searchedTerms]);

  // Genes for the Reported genes preset: the report's own, else the c.HGVS transcript's
  useEffect(() => {
    for (const record of reportedRecords) {
      const transcript = record.cHgvs?.split(':')[0];
      const lookupKey = `gene|${record.key}`;
      if (record.gene || !transcript || requestedLookupsRef.current.has(lookupKey)) {
        continue;
      }
      requestedLookupsRef.current.add(lookupKey);
      geneForTranscript(transcript).then((gene) => {
        if (gene) {
          setLookedUpGenes((previous) => new Map(previous).set(record.key, gene));
        }
      });
    }
  }, [reportedRecords]);
  // A fusion contributes each partner
  const reportedGenes = useMemo(() => [...new Set(reportedRecords
    .flatMap((record) => record.structural?.fusionGenes ?? [record.gene ?? lookedUpGenes.get(record.key)])
    .filter((gene): gene is string => Boolean(gene)))].sort(), [reportedRecords, lookedUpGenes]);

  const results = useMemo(
    () => layerReportedVariants(serverResults, searchedTerms, reportedRecords, normalizedVariants, fusionGeneRanges),
    [serverResults, searchedTerms, reportedRecords, normalizedVariants, fusionGeneRanges]
  );

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
    setServerResults([]);
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
    codeIndex: TumorTypeCodeIndex,
  ): Variant[] => {
    try {
      let filtered = applyFiltersToVariants(variants, filters);
      // Without a tumor type, "Actionable, this tumor type" stands as "Actionable, any tumor type"
      const actionableLabel = effectiveActionabilityFilter(filters.actionableFilter, Boolean(cancerType));

      if (actionableLabel === "Actionable, this tumor type") {
        // A-level evidence whose cancer-type code is allowed for the tumor type, or for all cancers
        filtered = filtered.map(variant => ({
          ...variant,
          txImplications: variant.txImplications?.filter(tx =>
            isALevelEvidence(tx.evidenceLevel) && matchesTumorType(tx.phenotypeCodings, cancerType, codeIndex)),
          dxImplications: variant.dxImplications?.filter(dx =>
            isALevelEvidence(dx.evidenceLevel) && matchesTumorType(dx.phenotypeCodings, cancerType, codeIndex))
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
      } else if (actionableLabel === "Reported variants") {
        // Only the patient's reported variants, with all their implications
        filtered = filtered.filter((variant) => (variant.reportedIn?.length ?? 0) > 0);
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
    const filtered = applyUnifiedFilters(results, filters, selectedCancerType, tumorTypeCodes);
    setFilteredResults(filtered);
  }, [results, selectedCancerType, tumorTypeCodes, applyUnifiedFilters]);

  // Apply filters whenever currentFilters, results, or related dependencies change
  useEffect(() => {
    if (results.length > 0) {
      const filtered = applyUnifiedFilters(results, currentFilters, selectedCancerType, tumorTypeCodes);
      setFilteredResults(filtered);
    } else {
      setFilteredResults([]);
    }
  }, [results, currentFilters, selectedCancerType, tumorTypeCodes, applyUnifiedFilters]);

  const hasActiveFilters = hasActiveSidebarFilters(currentFilters);
  const isSearching = Object.values(searchStatus).includes('searching');
  // Each variant records the search term it was found under
  const variantCountsByTerm = useMemo(() => {
    const counts: Record<string, number> = {};
    results.forEach((variant) => {
      counts[variant.range] = (counts[variant.range] ?? 0) + 1;
    });
    return counts;
  }, [results]);
  // A term the server found nothing in (or has no data for) still counts as found when a reported variant is there
  const displayedSearchStatus = useMemo(() => Object.fromEntries(Object.entries(searchStatus).map(([term, status]) =>
    [term, (status === 'no results' || status === 'not_on_server') && (variantCountsByTerm[term] ?? 0) > 0 ? 'found' : status])), [searchStatus, variantCountsByTerm]);
  // The genomics server has no data for this patient: results come only from the genetic report
  const noServerData = Object.values(searchStatus).includes('not_on_server');
  const collapsedSearchSummary = [
    searchInput.trim() || 'No search terms loaded',
    Object.keys(searchStatus).length > 0 ? searchSummaryText(displayedSearchStatus, noServerData) : undefined,
  ].filter(Boolean).join(' · ');

  const handleSearch = async () => {
    // Nothing to search until a patient is chosen
    if (!subjectId.trim()) {
      return;
    }
    if (!searchInput.trim()) {
      setServerResults([]);
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
      })
      // A term given twice ("ERBB2, erbb2") is searched once, so its variants aren't listed twice
      .filter((term, position, terms) => terms.indexOf(term) === position);
    // Cancel any ongoing searches
    abortControllersRef.current.forEach(controller => controller.abort());
    abortControllersRef.current.clear();

    // Create new abort controllers for this search
    const newControllers = new Map<string, AbortController>();
    const searchId = ++searchIdRef.current; // Use useRef for better ID generation

    setServerResults([]);
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

    // Once the server says it has no data for this patient, the remaining terms aren't sent
    let patientNotOnServer = false;

    // Process each search term sequentially and update results incrementally
    for (const term of searchTerms) {
      if (!term) continue;

      const controller = newControllers.get(term);
      if (!controller) continue;

      if (patientNotOnServer) {
        if (!controller.signal.aborted) {
          setSearchStatus(prev => ({ ...prev, [term]: 'not_on_server' }));
        }
        continue;
      }

      try {
        // Use the cached variant service - it handles all implications internally
        const processedVariants = await findSubjectVariantsWithCache(term, subjectId, searchId, enableCatVrsQueries);

        // Check if this search was cancelled
        if (controller.signal.aborted) {
          continue;
        }

        if (processedVariants.length > 0) {
          // Add all processed variants at once - they already have all implications
          setServerResults(prev => [...prev, ...processedVariants]);

          if (!controller.signal.aborted) {
            setSearchStatus(prev => ({ ...prev, [term]: 'found' }));
          }
        } else {
          if (!controller.signal.aborted) {
            setSearchStatus(prev => ({ ...prev, [term]: 'no results' }));
          }
        }
      } catch (error) {
        if (error instanceof PatientNotOnServerError) {
          // Not an error in the term: the patient's reported variants can still show for it
          patientNotOnServer = true;
          if (!controller.signal.aborted) {
            setSearchStatus(prev => ({ ...prev, [term]: 'not_on_server' }));
          }
        } else if (!controller.signal.aborted) {
          console.error(`Error searching for ${term}:`, error);
          setSearchStatus(prev => ({ ...prev, [term]: 'error' }));
          setInvalidRanges(prev => [...prev, term]);
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
          {CASE_VIEWS.map(({ view, label, icon: Icon, iconClass, activeIconClass, tooltip }) => {
            const isActive = caseView === view;
            return (
              <button
                key={view}
                type="button"
                role="tab"
                title={tooltip}
                aria-selected={isActive}
                onClick={() => chooseCaseView(view)}
                className={`inline-flex items-center gap-2 rounded-lg border px-5 py-2 text-sm font-semibold ${isActive ? 'border-slate-900 bg-slate-900 text-white' : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'}`}
              >
                <Icon size={16} aria-hidden="true" className={isActive ? activeIconClass : iconClass} />
                {label}
              </button>
            );
          })}
        </div>

        {caseView === 'summary' && (
          // White cards on a light gray ground give the summary some depth
          <div className="mb-8 rounded-3xl bg-slate-100 p-3 ring-1 ring-inset ring-slate-200 sm:p-5">
            <CaseSummary
              onGoToGenomicFindings={() => chooseCaseView('genomic')}
              subjectId={subjectId}
              clinicalState={clinicalState}
              onRetry={retryClinicalRecord}
              crosswalk={tumorTypeCrosswalk}
            />
          </div>
        )}

        {/* Placeholder until imaging data (slides, diagnostic images) is pulled in */}
        {caseView === 'imaging' && (
          <div className="mb-8 rounded-3xl bg-slate-100 p-3 ring-1 ring-inset ring-slate-200 sm:p-5">
            <div className="rounded-2xl border border-slate-200 bg-white p-8 text-center text-slate-600 shadow-sm">
              No images available.
            </div>
          </div>
        )}

        {/* Kept mounted while hidden, so a generated synthesis is still there when you return to the tab */}
        <div className={caseView === 'synthesis' ? 'mb-8 rounded-3xl bg-slate-100 p-3 ring-1 ring-inset ring-slate-200 sm:p-5' : 'hidden'}>
          <SynthesisView subjectId={subjectId} />
        </div>

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
          {/* A light emerald wash (the Genomics Explorer tab's color) makes the search box stand out on the white page */}
          <section aria-label="Variant search" className="mb-8 rounded-2xl border border-emerald-200 bg-gradient-to-b from-emerald-50 to-white p-5 shadow-[0_1px_3px_rgba(15,23,42,0.06),0_4px_12px_rgba(15,23,42,0.04)]">
            <div className={`flex flex-wrap items-center justify-between gap-3 ${isSearchWorkspaceCollapsed ? '' : 'mb-4 border-b border-emerald-100 pb-3'}`}>
              <div className="flex min-w-0 flex-1 items-center gap-2.5">
                <span className="flex h-8 w-8 shrink-0 items-center justify-center rounded-lg bg-emerald-100 text-emerald-700" aria-hidden="true">
                  <Dna size={16} strokeWidth={2.25} />
                </span>
                <h2 className="shrink-0 text-base font-semibold tracking-tight text-slate-900">Variant search</h2>
                {/* Hidden, the box shrinks to one line: the search terms and how the search went */}
                {isSearchWorkspaceCollapsed && (
                  <span className="min-w-0 flex-1 truncate text-sm text-slate-600" title={collapsedSearchSummary}>{collapsedSearchSummary}</span>
                )}
              </div>
              <button
                type="button"
                onClick={() => setIsSearchWorkspaceCollapsed((currentValue) => !currentValue)}
                aria-expanded={!isSearchWorkspaceCollapsed}
                className="inline-flex items-center gap-1 text-sm font-medium text-slate-600 hover:text-slate-900"
              >
                {isSearchWorkspaceCollapsed ? <><ChevronDown size={16} aria-hidden="true" /> Show</> : <><ChevronUp size={16} aria-hidden="true" /> Hide</>}
              </button>
            </div>

            {!isSearchWorkspaceCollapsed && (
              <div className="flex flex-col gap-3">
                {/* Sets this tab apart from the lab's report on the Case Summary: this is the complete data set */}
                <p className="-mt-1 text-sm text-slate-600">
                  Search the patient&apos;s complete genomic data set, including variants not in the curated report.
                </p>
                <SearchForm
                  searchInput={searchInput}
                  setSearchInput={handleSearchInputChange}
                  handleSearch={handleSearch}
                  isSearching={isSearching}
                  selectedCancerType={selectedCancerType}
                  onPresetSelect={handlePresetSelect}
                  reportedGenes={reportedGenes}
                  onReportedGenesSelect={() => handleSearchInputChange(reportedGenes.join(', '))}
                  enableCatVrsQueries={enableCatVrsQueries}
                  onEnableCatVrsQueriesChange={setEnableCatVrsQueries}
                />
                <SearchStatus
                  searchStatus={displayedSearchStatus}
                  noServerData={noServerData}
                  invalidRanges={invalidRanges}
                  variantCounts={variantCountsByTerm}
                  onCancelSearch={handleCancelSearch}
                />
              </div>
            )}
          </section>
          <FilterSidebar
            onFilterChange={handleFilterChange}
            isOpen={isFilterOpen}
            onClose={() => setIsFilterOpen(false)}
            results={results}
            currentFilters={currentFilters}
            hasTumorType={Boolean(selectedCancerType)}
          />
          {/* Folder tabs: the active one is white with the emerald outline and joins the card below; the other has a
              light outline and fill so it reads as clickable */}
          <div ref={tabBarRef} className="relative z-10 -mb-px flex flex-wrap gap-1.5 pl-4" role="tablist">
            {(['results', 'pathways'] as const).map((view) => {
              const isActive = activeView === view;
              const isDisabled = view === 'pathways' && selectedVariants.length === 0;
              const TabIcon = view === 'results' ? Table2 : Network;
              const count = view === 'results' ? String(filteredResults.length) : selectedVariants.length ? `${selectedVariants.length} selected` : undefined;
              return (
                <button
                  key={view}
                  type="button"
                  role="tab"
                  aria-selected={isActive}
                  disabled={isDisabled}
                  title={isDisabled ? 'Select variants in the results table first' : undefined}
                  onClick={() => switchView(view)}
                  className={`inline-flex items-center gap-2 rounded-t-xl border px-4 py-2 text-sm ${isActive ? 'border-emerald-200 border-b-white bg-white font-semibold text-slate-900' : 'border-slate-200 border-b-emerald-200 bg-slate-50 font-medium text-slate-600 hover:bg-slate-100'} disabled:cursor-not-allowed disabled:opacity-50`}
                >
                  <TabIcon size={15} aria-hidden="true" className={isActive ? 'text-emerald-700' : 'text-slate-500'} />
                  {view === 'results' ? 'Results' : 'Pathways'}
                  {count && (
                    <span className={`rounded-full px-2 py-0.5 text-xs ${isActive ? 'bg-emerald-50 text-emerald-800' : 'bg-slate-200 text-slate-600'}`}>{count}</span>
                  )}
                </button>
              );
            })}
          </div>
          {/* The results table stays mounted so its state (e.g., oncogenicity results) survives switching tabs */}
          <div className={activeView === 'results' ? '' : 'hidden'}>
            <ResultsTable
              results={filteredResults}
              hasSearched={Object.keys(searchStatus).length > 0}
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
            <div className="min-h-screen rounded-2xl border border-emerald-200 bg-white p-3 shadow-[0_1px_3px_rgba(15,23,42,0.06),0_4px_12px_rgba(15,23,42,0.04)] sm:p-6">
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
