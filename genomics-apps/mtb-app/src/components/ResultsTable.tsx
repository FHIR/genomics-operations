"use client";

import { MouseEvent as ReactMouseEvent, useEffect, useMemo, useRef, useState } from 'react';
import Image from 'next/image';
import { ArrowDown, ArrowDownToLine, ArrowUp, ArrowUpDown, Dna, Filter, Search, SlidersHorizontal } from 'lucide-react';
import { Variant } from '@/types/variants';
import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import DxImplicationCell from './DxImplicationCell';
import MolecularConsequenceCell from './MolecularConsequenceCell';
import TxImplicationCell from './TxImplicationCell';
import { CooccurrenceProfile } from '@/services/cooccurrenceService';
import VariantGroupRow from './VariantGroupRow';
import {
  ColumnId,
  DEFAULT_COLUMN_ORDER,
  RESULTS_TABLE_COLUMNS,
} from './resultsTableColumns';
import {
  getConceptDisplayText,
  getConceptNarrativeText,
  getGaugeImagePath,
  getGaugeKindForScore,
  getObservationCaveat,
  ONCOGENICITY_RETRY_DELAYS_MS,
  ONCOGENICITY_SLOW_REQUEST_MS,
  OncogenicityPredictionResult,
  runWithRetries,
} from '@/types/oncogenicity';
import { useHasElapsed } from '@/lib/useHasElapsed';
import _ from 'lodash';

type SortDirection = 'asc' | 'desc';

interface SortState {
  columnId: ColumnId;
  direction: SortDirection;
}

const DEFAULT_SORT_STATE: SortState = {
  columnId: 'range',
  direction: 'asc',
};

const COLUMN_MAP = Object.fromEntries(
  RESULTS_TABLE_COLUMNS.map((column) => [column.id, column])
) as Record<ColumnId, (typeof RESULTS_TABLE_COLUMNS)[number]>;

const DEFAULT_COLUMN_VISIBILITY = Object.fromEntries(
  RESULTS_TABLE_COLUMNS.map((column) => [column.id, !['genomicSourceClass', 'variantAlleleFrequency', 'oncogenicityPrediction'].includes(column.id)])
) as Record<ColumnId, boolean>;

const DEFAULT_COLUMN_WIDTHS = Object.fromEntries(
  RESULTS_TABLE_COLUMNS.map((column) => [column.id, column.defaultWidth])
) as Record<ColumnId, number>;

const TABLE_SETTINGS_STORAGE_KEY = 'mtb-results-table-settings';
const TABLE_WIDTHS_VERSION = 6;
const TABLE_COLUMN_VISIBILITY_VERSION = 2;
const DEFAULT_ENABLE_CAT_VRS_QUERIES = false;

function normalizeColumnOrder(columnOrder?: ColumnId[]) {
  const knownColumnIds = new Set(DEFAULT_COLUMN_ORDER);
  const validOrder = (columnOrder ?? []).filter(
    (columnId): columnId is ColumnId => knownColumnIds.has(columnId)
  );
  const normalizedOrder = [...new Set(validOrder)];

  for (const defaultColumnId of DEFAULT_COLUMN_ORDER) {
    if (normalizedOrder.includes(defaultColumnId)) {
      continue;
    }

    const nextExistingDefaultColumn = DEFAULT_COLUMN_ORDER
      .slice(DEFAULT_COLUMN_ORDER.indexOf(defaultColumnId) + 1)
      .find((columnId) => normalizedOrder.includes(columnId));

    if (!nextExistingDefaultColumn) {
      normalizedOrder.push(defaultColumnId);
      continue;
    }

    normalizedOrder.splice(normalizedOrder.indexOf(nextExistingDefaultColumn), 0, defaultColumnId);
  }

  return normalizedOrder;
}

interface ResultsTableProps {
  results: Variant[];
  // Whether a search has run for this patient, for the empty-state message
  hasSearched: boolean;
  predictorTumorType?: string;
  onToggleFilters: () => void;
  hasActiveFilters: boolean;
  enableCatVrsQueries: boolean;
  onEnableCatVrsQueriesChange: (enabled: boolean) => void;
  cooccurrencesByVariant: Record<string, CooccurrenceProfile[]>;
  selectedVariantIds: Set<string>;
  onToggleVariantSelection: (variantId: string) => void;
  onSelectVariants: (variantIds: string[]) => void;
  onVisibleVariantIdsChange: (variantIds: string[]) => void;
  highlightRequest: HighlightRequest | null;
  onHighlightMissing: (variantId: string) => void;
}

export interface HighlightRequest {
  variantId: string;
  requestId: number;
}

const SELECTION_COLUMN_WIDTH = 44;

function getSimpleVariantLabel(variantString: string) {
  const parts = variantString.split(':');

  if (parts.length !== 4) {
    return 'Simple';
  }

  const [, , deleted = '', inserted = ''] = parts;

  if (deleted.length === inserted.length) {
    if (deleted.length === 1) {
      return 'SNV';
    }

    if (deleted.length > 1) {
      return 'MNV';
    }
  }

  if (deleted.length !== inserted.length) {
    return 'InDel';
  }

  return 'Simple';
}

function getVariantSortValue(variant: Variant) {
  if (variant.variantType === 'simple') {
    return `${getSimpleVariantLabel(variant.variant)} ${variant.variant}`;
  }

  return variant.variant;
}

async function getResponseErrorMessage(response: Response, fallbackMessage: string) {
  try {
    const payload = await response.json() as { error?: string };
    return payload.error || `${fallbackMessage} (${response.status})`;
  } catch {
    return `${fallbackMessage} (${response.status})`;
  }
}

// Flat outline buttons for the results toolbar; the dark primary style is kept for Search
const TOOLBAR_BUTTON_CLASSES = 'inline-flex items-center gap-2 rounded-lg border border-slate-300 bg-white px-3.5 py-2 text-sm font-semibold text-slate-700 transition-colors hover:bg-slate-50 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2';

function getOncogenicityKey(variant: Variant) {
  return variant.id ?? variant.variant;
}

function getColumnTextValue(
  variant: Variant,
  columnId: ColumnId,
  range: string,
  oncogenicityResults?: Record<string, OncogenicityPredictionResult>,
) {
  switch (columnId) {
    case 'range':
      return range;
    case 'variant':
      return [getVariantSortValue(variant), variant.molecularConsequences?.[0]?.proteinChange]
        .filter(Boolean)
        .join(' ');
    case 'genomicSourceClass':
      return variant.genomicSourceClass ?? '<unknown>';
    case 'variantAlleleFrequency':
      return variant.variantAlleleFrequency?.toString() ?? '<none found>';
    case 'oncogenicityPrediction':
      return [
        oncogenicityResults?.[getOncogenicityKey(variant)]?.score,
        oncogenicityResults?.[getOncogenicityKey(variant)]?.interpretation,
        variant.oncogenicityPrediction,
      ]
        .filter((value) => value !== undefined && value !== null && String(value).trim() !== '')
        .join(' ');
    case 'molecularConsequences':
      return (variant.molecularConsequences ?? [])
        .flatMap((consequence) => [consequence.impact, consequence.featureConsequence, consequence.proteinChange])
        .filter(Boolean)
        .join(' ');
    case 'dxImplications':
      return (variant.dxImplications ?? [])
        .flatMap((implication) => [
          implication.clinicalSignificance,
          implication.predictedPhenotype,
          implication.evidenceLevel,
          implication.variantId,
        ])
        .filter(Boolean)
        .join(' ');
    case 'txImplications':
      return (variant.txImplications ?? [])
        .flatMap((implication) => [
          implication.therapeuticImplication,
          implication.medication,
          implication.evidenceLevel,
          implication.phenotypicContext,
        ])
        .filter(Boolean)
        .join(' ');
  }
}

function compareColumnValues(left: string, right: string, direction: SortDirection) {
  const normalizedLeft = left.trim().toLocaleLowerCase();
  const normalizedRight = right.trim().toLocaleLowerCase();

  if (!normalizedLeft && !normalizedRight) {
    return 0;
  }

  if (!normalizedLeft) {
    return 1;
  }

  if (!normalizedRight) {
    return -1;
  }

  const leftNumber = Number(normalizedLeft);
  const rightNumber = Number(normalizedRight);

  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) {
    return direction === 'asc' ? leftNumber - rightNumber : rightNumber - leftNumber;
  }

  const comparison = normalizedLeft.localeCompare(normalizedRight, undefined, {
    numeric: true,
    sensitivity: 'base',
  });

  return direction === 'asc' ? comparison : -comparison;
}

function getMinimumColumnWidth(columnId: ColumnId) {
  return COLUMN_MAP[columnId]?.minWidth ?? 80;
}

export default function ResultsTable({
  results,
  hasSearched,
  predictorTumorType,
  onToggleFilters,
  hasActiveFilters,
  enableCatVrsQueries,
  onEnableCatVrsQueriesChange,
  cooccurrencesByVariant,
  selectedVariantIds,
  onToggleVariantSelection,
  onSelectVariants,
  onVisibleVariantIdsChange,
  highlightRequest,
  onHighlightMissing,
}: ResultsTableProps) {
  const [highlightedVariantId, setHighlightedVariantId] = useState<string | null>(null);
  const [columnOrder, setColumnOrder] = useState<ColumnId[]>(DEFAULT_COLUMN_ORDER);
  const [columnVisibility, setColumnVisibility] = useState<Record<ColumnId, boolean>>(
    DEFAULT_COLUMN_VISIBILITY
  );
  const [columnWidths, setColumnWidths] = useState<Record<ColumnId, number>>(
    DEFAULT_COLUMN_WIDTHS
  );
  const [draggedColumnId, setDraggedColumnId] = useState<ColumnId | null>(null);
  const [sortState, setSortState] = useState<SortState | null>(DEFAULT_SORT_STATE);
  const [isCustomizeTableOpen, setIsCustomizeTableOpen] = useState(false);
  const [oncogenicityResults, setOncogenicityResults] = useState<Record<string, OncogenicityPredictionResult>>({});
  const [selectedOncogenicityKey, setSelectedOncogenicityKey] = useState<string | null>(null);
  const [showClinVarWorkInProgress, setShowClinVarWorkInProgress] = useState(false);
  const [columnFilters, setColumnFilters] = useState<Record<ColumnId, string>>({
    range: '',
    variant: '',
    genomicSourceClass: '',
    variantAlleleFrequency: '',
    oncogenicityPrediction: '',
    molecularConsequences: '',
    dxImplications: '',
    txImplications: '',
  });
  const customizeTableRef = useRef<HTMLDivElement | null>(null);
  const extendedEvidenceRef = useRef<HTMLElement | null>(null);
  const pendingEvidenceScrollKeyRef = useRef<string | null>(null);
  // Predictions depend on the tumor type, so a Dx change discards them, including any still running
  const tumorTypeGenerationRef = useRef(0);

  useEffect(() => {
    tumorTypeGenerationRef.current += 1;
    pendingEvidenceScrollKeyRef.current = null;
    setOncogenicityResults({});
    setSelectedOncogenicityKey(null);
  }, [predictorTumorType]);

  useEffect(() => {
    if (typeof window === 'undefined') {
      return;
    }

    const storedSettings = window.localStorage.getItem(TABLE_SETTINGS_STORAGE_KEY);

    if (!storedSettings) {
      return;
    }

    try {
      const parsedSettings = JSON.parse(storedSettings) as {
        columnOrder?: ColumnId[];
        columnVisibility?: Partial<Record<ColumnId, boolean>>;
        columnVisibilityVersion?: number;
        columnWidths?: Partial<Record<ColumnId, number>>;
        columnWidthsVersion?: number;
        enableCatVrsQueries?: boolean;
      };

      setColumnOrder(normalizeColumnOrder(parsedSettings.columnOrder));
      const persistedVisibility = parsedSettings.columnVisibilityVersion === TABLE_COLUMN_VISIBILITY_VERSION
        ? parsedSettings.columnVisibility
        : {};

      setColumnVisibility({
        ...DEFAULT_COLUMN_VISIBILITY,
        ...persistedVisibility,
      });
      const persistedWidths = parsedSettings.columnWidthsVersion === TABLE_WIDTHS_VERSION
        ? Object.fromEntries(
          Object.entries(parsedSettings.columnWidths ?? {}).map(([columnId, width]) => {
            const typedColumnId = columnId as ColumnId;
            const minimumWidth = getMinimumColumnWidth(typedColumnId);

            return [typedColumnId, Math.max(Number(width) || minimumWidth, minimumWidth)];
          })
        )
        : {};

      setColumnWidths({
        ...DEFAULT_COLUMN_WIDTHS,
        ...persistedWidths,
      });
      onEnableCatVrsQueriesChange(parsedSettings.enableCatVrsQueries ?? DEFAULT_ENABLE_CAT_VRS_QUERIES);
    } catch {
      window.localStorage.removeItem(TABLE_SETTINGS_STORAGE_KEY);
    }
  }, [onEnableCatVrsQueriesChange]);

  useEffect(() => {
    if (typeof window === 'undefined') {
      return;
    }

    window.localStorage.setItem(
      TABLE_SETTINGS_STORAGE_KEY,
      JSON.stringify({
        columnOrder,
        columnVisibility,
        columnVisibilityVersion: TABLE_COLUMN_VISIBILITY_VERSION,
        columnWidths,
        columnWidthsVersion: TABLE_WIDTHS_VERSION,
        enableCatVrsQueries,
      })
    );
  }, [columnOrder, columnVisibility, columnWidths, enableCatVrsQueries]);

  useEffect(() => {
    if (!isCustomizeTableOpen) {
      return;
    }

    const handlePointerDown = (event: MouseEvent) => {
      if (!customizeTableRef.current?.contains(event.target as Node)) {
        setIsCustomizeTableOpen(false);
      }
    };

    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        setIsCustomizeTableOpen(false);
      }
    };

    document.addEventListener('mousedown', handlePointerDown);
    document.addEventListener('keydown', handleKeyDown);

    return () => {
      document.removeEventListener('mousedown', handlePointerDown);
      document.removeEventListener('keydown', handleKeyDown);
    };
  }, [isCustomizeTableOpen]);

  // Function to render Dx Implications content
  const renderDxImplications = (implications?: DxImplication[]) => {
    return <DxImplicationCell implications={implications} />;
  };

  // Function to render Molecular Consequences content
  const renderMolecularConsequences = (consequences?: MolecularConsequence[]) => {
    return <MolecularConsequenceCell consequences={consequences} />;
  };

  // Function to render Tx Implications content
  const renderTxImplications = (implications: ProcessedTxImplication[] | undefined, variant: Variant) => {
    return (
      <TxImplicationCell
        implications={implications}
        cooccurrences={variant.id ? cooccurrencesByVariant[variant.id] : undefined}
      />
    );
  };

  const visibleColumns = useMemo(
    () => columnOrder.filter((columnId) => columnVisibility[columnId]).map((columnId) => COLUMN_MAP[columnId]),
    [columnOrder, columnVisibility]
  );

  const groupedResults = useMemo(
    () => Object.entries(_.groupBy(results, 'range')),
    [results]
  );

  const processedGroups = useMemo(() => {
    const activeSortState = sortState ?? DEFAULT_SORT_STATE;
    const activeVisibleFilters = visibleColumns.filter(
      (column) => columnFilters[column.id].trim() !== ''
    );

    const filteredGroups = groupedResults
      .map(([range, variants]) => {
        const filteredVariants = variants.filter((variant) => {
          return activeVisibleFilters.every((column) => {
            const filterValue = columnFilters[column.id].trim().toLocaleLowerCase();
            const cellValue = getColumnTextValue(variant, column.id, range, oncogenicityResults).toLocaleLowerCase();
            return cellValue.includes(filterValue);
          });
        });

        return [range, filteredVariants] as const;
      })
      .filter(([, variants]) => variants.length > 0);

    const sortedGroups = filteredGroups
      .map(([range, variants]) => {
        const sortedVariants = [...variants].sort((leftVariant, rightVariant) => {
          const leftValue = getColumnTextValue(leftVariant, activeSortState.columnId, range, oncogenicityResults);
          const rightValue = getColumnTextValue(rightVariant, activeSortState.columnId, range, oncogenicityResults);
          return compareColumnValues(leftValue, rightValue, activeSortState.direction);
        });

        return [range, sortedVariants] as const;
      })
      .sort(([leftRange, leftVariants], [rightRange, rightVariants]) => {
        const leftValue = getColumnTextValue(leftVariants[0], activeSortState.columnId, leftRange, oncogenicityResults);
        const rightValue = getColumnTextValue(rightVariants[0], activeSortState.columnId, rightRange, oncogenicityResults);
        return compareColumnValues(leftValue, rightValue, activeSortState.direction);
      });

    return sortedGroups;
  }, [groupedResults, visibleColumns, columnFilters, sortState, oncogenicityResults]);

  const visibleVariantIds = useMemo(
    () => processedGroups.flatMap(([, variants]) => variants.map((variant) => variant.id).filter((id): id is string => Boolean(id))),
    [processedGroups]
  );

  useEffect(() => {
    onVisibleVariantIdsChange(visibleVariantIds);
  }, [visibleVariantIds, onVisibleVariantIdsChange]);

  // "Show in results table": scroll to the row and highlight it briefly
  useEffect(() => {
    if (!highlightRequest) {
      return;
    }

    const { variantId } = highlightRequest;

    if (!visibleVariantIds.includes(variantId)) {
      onHighlightMissing(variantId);
      return;
    }

    setHighlightedVariantId(variantId);
    const scrollTimer = window.setTimeout(() => {
      document.getElementById(`variant-row-${variantId}`)?.scrollIntoView({ behavior: 'smooth', block: 'center' });
    }, 80);
    const clearTimer = window.setTimeout(() => setHighlightedVariantId(null), 2400);

    return () => {
      window.clearTimeout(scrollTimer);
      window.clearTimeout(clearTimer);
    };
    // Only react to new requests, not to later changes in visible rows
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [highlightRequest]);

  const totalTableWidth = SELECTION_COLUMN_WIDTH + visibleColumns.reduce(
    (width, column) => width + (columnWidths[column.id] ?? column.defaultWidth),
    0
  );
  const txColumnWidth = columnWidths.txImplications ?? COLUMN_MAP.txImplications.defaultWidth;
  const hasExplicitTxColumnWidth = txColumnWidth !== COLUMN_MAP.txImplications.defaultWidth;
  const minimumTableWidth = SELECTION_COLUMN_WIDTH + visibleColumns.reduce(
    (width, column) => width + getMinimumColumnWidth(column.id),
    0
  );

  const hasActiveTableFilters = visibleColumns.some((column) => columnFilters[column.id].trim() !== '');

  const handleColumnVisibilityChange = (columnId: ColumnId) => {
    setColumnVisibility((currentVisibility) => {
      const visibleColumnCount = Object.values(currentVisibility).filter(Boolean).length;

      if (currentVisibility[columnId] && visibleColumnCount === 1) {
        return currentVisibility;
      }

      return {
        ...currentVisibility,
        [columnId]: !currentVisibility[columnId],
      };
    });
  };

  const moveColumn = (columnId: ColumnId, direction: 'left' | 'right') => {
    setColumnOrder((currentOrder) => {
      const columnIndex = currentOrder.indexOf(columnId);
      const targetIndex = direction === 'left' ? columnIndex - 1 : columnIndex + 1;

      if (columnIndex === -1 || targetIndex < 0 || targetIndex >= currentOrder.length) {
        return currentOrder;
      }

      const updatedOrder = [...currentOrder];
      const [column] = updatedOrder.splice(columnIndex, 1);
      updatedOrder.splice(targetIndex, 0, column);
      return updatedOrder;
    });
  };

  const resetColumns = () => {
    setColumnOrder(DEFAULT_COLUMN_ORDER);
    setColumnVisibility(DEFAULT_COLUMN_VISIBILITY);
    setColumnWidths(DEFAULT_COLUMN_WIDTHS);
    onEnableCatVrsQueriesChange(DEFAULT_ENABLE_CAT_VRS_QUERIES);
    setSortState(DEFAULT_SORT_STATE);
    setColumnFilters({
      range: '',
      variant: '',
      genomicSourceClass: '',
      variantAlleleFrequency: '',
      oncogenicityPrediction: '',
      molecularConsequences: '',
      dxImplications: '',
      txImplications: '',
    });
  };

  const toggleSort = (columnId: ColumnId) => {
    setSortState((currentSortState) => {
      if (!currentSortState || currentSortState.columnId !== columnId) {
        return { columnId, direction: 'asc' };
      }

      if (currentSortState.direction === 'asc') {
        return { columnId, direction: 'desc' };
      }

      return null;
    });
  };

  const clearColumnFilters = () => {
    setColumnFilters((currentFilters) => Object.fromEntries(Object.keys(currentFilters).map((columnId) => [columnId, ''])) as Record<ColumnId, string>);
  };

  const updateColumnFilter = (columnId: ColumnId, value: string) => {
    setColumnFilters((currentFilters) => ({
      ...currentFilters,
      [columnId]: value,
    }));
  };

  const handleDrop = (targetColumnId: ColumnId) => {
    if (!draggedColumnId || draggedColumnId === targetColumnId) {
      setDraggedColumnId(null);
      return;
    }

    setColumnOrder((currentOrder) => {
      const draggedIndex = currentOrder.indexOf(draggedColumnId);
      const targetIndex = currentOrder.indexOf(targetColumnId);

      if (draggedIndex === -1 || targetIndex === -1) {
        return currentOrder;
      }

      const updatedOrder = [...currentOrder];
      updatedOrder.splice(draggedIndex, 1);
      updatedOrder.splice(targetIndex, 0, draggedColumnId);
      return updatedOrder;
    });
    setDraggedColumnId(null);
  };

  const startResize = (event: ReactMouseEvent<HTMLDivElement>, columnId: ColumnId) => {
    event.preventDefault();
    event.stopPropagation();

    const startingX = event.clientX;
    const startingWidth = columnWidths[columnId];
    const minimumWidth = getMinimumColumnWidth(columnId);

    const handleMouseMove = (mouseEvent: MouseEvent) => {
      const nextWidth = Math.max(startingWidth + (mouseEvent.clientX - startingX), minimumWidth);

      setColumnWidths((currentWidths) => ({
        ...currentWidths,
        [columnId]: nextWidth,
      }));
    };

    const handleMouseUp = () => {
      window.removeEventListener('mousemove', handleMouseMove);
      window.removeEventListener('mouseup', handleMouseUp);
    };

    window.addEventListener('mousemove', handleMouseMove);
    window.addEventListener('mouseup', handleMouseUp);
  };

  const selectedOncogenicityResult = selectedOncogenicityKey ? oncogenicityResults[selectedOncogenicityKey] : undefined;
  const isEvidenceSlow = useHasElapsed(
    selectedOncogenicityResult?.evidenceStatus === 'loading' ? selectedOncogenicityResult.evidenceStartedAt : undefined,
    ONCOGENICITY_SLOW_REQUEST_MS,
  );

  // Scroll to the extended evidence once it arrives for the request the user started
  useEffect(() => {
    if (
      selectedOncogenicityResult?.evidenceStatus === 'ready'
      && pendingEvidenceScrollKeyRef.current === selectedOncogenicityResult.key
    ) {
      pendingEvidenceScrollKeyRef.current = null;
      extendedEvidenceRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  }, [selectedOncogenicityResult?.evidenceStatus, selectedOncogenicityResult?.key]);

  const computeOncogenicityPrediction = async (variant: Variant) => {
    const generation = tumorTypeGenerationRef.current;
    const setResultsIfCurrent: typeof setOncogenicityResults = (update) => {
      if (tumorTypeGenerationRef.current === generation) {
        setOncogenicityResults(update);
      }
    };
    const key = getOncogenicityKey(variant);

    if (oncogenicityResults[key]?.status === 'loading') {
      return;
    }

    setResultsIfCurrent((currentResults) => ({
      ...currentResults,
      [key]: {
        key,
        variant: variant.variant,
        status: 'loading',
        startedAt: Date.now(),
        gauge: 'undetermined',
        evidenceStatus: 'idle',
      },
    }));

    const requestPrediction = async () => {
      const response = await fetch('/api/oncogenicity/predict', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          variant: variant.variant,
          tumorType: predictorTumorType,
        }),
      });

      if (!response.ok) {
        throw new Error(await getResponseErrorMessage(response, 'Prediction request failed'));
      }

      const payload = await response.json() as { observation?: OncogenicityPredictionResult['observation'] };
      const rawScore = payload.observation?.valueInteger;
      const normalizedScore = typeof rawScore === 'number' ? rawScore : Number(rawScore);
      const score = Number.isFinite(normalizedScore) ? normalizedScore : undefined;

      // A response without a score is treated as a failure; the predictor's annotation step fails intermittently
      if (score === undefined) {
        throw new Error(
          getConceptNarrativeText(payload.observation?.dataAbsentReason) ?? 'Prediction did not return a score'
        );
      }

      return { payload, score };
    };

    try {
      const { payload, score } = await runWithRetries(requestPrediction, (retryAttempt, error) => {
        setResultsIfCurrent((currentResults) => ({
          ...currentResults,
          [key]: {
            ...currentResults[key],
            retryAttempt,
            errorMessage: error instanceof Error ? error.message : undefined,
          },
        }));
      });

      setResultsIfCurrent((currentResults) => ({
        ...currentResults,
        [key]: {
          key,
          variant: variant.variant,
          status: 'ready',
          gauge: getGaugeKindForScore(score),
          score,
          interpretation: getConceptDisplayText(payload.observation?.interpretation),
          observation: payload.observation,
          evidenceStatus: 'idle',
        },
      }));
    } catch (error) {
      setResultsIfCurrent((currentResults) => ({
        ...currentResults,
        [key]: {
          key,
          variant: variant.variant,
          status: 'error',
          gauge: 'undetermined',
          evidenceStatus: 'idle',
          errorMessage: error instanceof Error ? error.message : 'Prediction unavailable',
        },
      }));
    }
  };

  const openOncogenicityDetails = (variant: Variant) => {
    setSelectedOncogenicityKey(getOncogenicityKey(variant));
    setShowClinVarWorkInProgress(false);
  };

  const scrollToExtendedEvidence = () => {
    extendedEvidenceRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
  };

  const loadExtendedEvidence = async () => {
    const generation = tumorTypeGenerationRef.current;
    const setResultsIfCurrent: typeof setOncogenicityResults = (update) => {
      if (tumorTypeGenerationRef.current === generation) {
        setOncogenicityResults(update);
      }
    };
    if (!selectedOncogenicityResult || selectedOncogenicityResult.evidenceStatus === 'loading') {
      return;
    }

    if (selectedOncogenicityResult.evidenceStatus === 'ready') {
      scrollToExtendedEvidence();
      return;
    }

    const key = selectedOncogenicityResult.key;
    pendingEvidenceScrollKeyRef.current = key;

    setResultsIfCurrent((currentResults) => ({
      ...currentResults,
      [key]: {
        ...currentResults[key],
        evidenceStatus: 'loading',
        evidenceStartedAt: Date.now(),
        evidenceRetryAttempt: undefined,
        evidenceError: undefined,
      },
    }));

    const requestEvidence = async () => {
      const response = await fetch('/api/oncogenicity/evidence', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          variant: selectedOncogenicityResult.variant,
          tumorType: predictorTumorType,
        }),
      });

      if (!response.ok) {
        throw new Error(await getResponseErrorMessage(response, 'Extended evidence request failed'));
      }

      return response.json() as Promise<{ evidence?: unknown }>;
    };

    try {
      const payload = await runWithRetries(requestEvidence, (retryAttempt) => {
        setResultsIfCurrent((currentResults) => ({
          ...currentResults,
          [key]: {
            ...currentResults[key],
            evidenceRetryAttempt: retryAttempt,
          },
        }));
      });

      setResultsIfCurrent((currentResults) => ({
        ...currentResults,
        [key]: {
          ...currentResults[key],
          evidenceStatus: 'ready',
          evidenceRetryAttempt: undefined,
          evidenceJson: payload.evidence,
          evidenceError: undefined,
        },
      }));
    } catch (error) {
      setResultsIfCurrent((currentResults) => ({
        ...currentResults,
        [key]: {
          ...currentResults[key],
          evidenceStatus: 'error',
          evidenceRetryAttempt: undefined,
          evidenceError: error instanceof Error ? error.message : 'Unable to load extended evidence',
        },
      }));
    }
  };

  const evidenceStatus = selectedOncogenicityResult?.evidenceStatus;
  const evidenceButtonLabel = evidenceStatus === 'loading'
    ? selectedOncogenicityResult?.evidenceRetryAttempt
      ? `Retrying (${selectedOncogenicityResult.evidenceRetryAttempt} of ${ONCOGENICITY_RETRY_DELAYS_MS.length})...`
      : isEvidenceSlow
        ? 'Waking up prediction server, this can take up to a minute...'
        : 'Loading evidence...'
    : evidenceStatus === 'ready'
      ? 'Evidence loaded below'
      : evidenceStatus === 'error'
        ? 'Retry loading evidence'
        : 'View extended evidence details';

  return (
    // One white card with the Genomic Findings tab's emerald as its border; joins the Results tab above it
    <div className="rounded-2xl border border-emerald-200 bg-white p-3 shadow-[0_1px_3px_rgba(15,23,42,0.06),0_4px_12px_rgba(15,23,42,0.04)] sm:p-4">
      <div className="mb-3 flex flex-wrap items-center justify-between gap-3">
        <div className="flex flex-wrap items-center gap-2">
          <button
            type="button"
            onClick={onToggleFilters}
            className={TOOLBAR_BUTTON_CLASSES}
          >
            <Filter size={15} aria-hidden="true" />
            Filters
            {hasActiveFilters && (
              <span className="rounded-full bg-amber-100 px-2 py-0.5 text-xs font-semibold text-amber-800 ring-1 ring-inset ring-amber-200">
                Active
              </span>
            )}
          </button>
          <div className="relative" ref={customizeTableRef}>
            <button
              type="button"
              aria-expanded={isCustomizeTableOpen}
              aria-haspopup="dialog"
              onClick={() => setIsCustomizeTableOpen((currentValue) => !currentValue)}
              className={TOOLBAR_BUTTON_CLASSES}
            >
              <SlidersHorizontal size={15} aria-hidden="true" />
              Customize table
            </button>
            {isCustomizeTableOpen && (
              <div className="absolute left-0 z-30 mt-2 w-96 rounded-lg border border-gray-200 bg-white p-4 shadow-xl">
                <div className="space-y-3">
                  {RESULTS_TABLE_COLUMNS.map((column) => {
                    const columnId = column.id;
                    const index = columnOrder.indexOf(columnId);
                    const isOnlyVisibleColumn = columnVisibility[columnId] && visibleColumns.length === 1;

                    return (
                      <div key={columnId} className="flex items-center justify-between gap-3 rounded-md border border-gray-100 px-3 py-2">
                        <label className="flex items-center gap-2 text-sm text-gray-700">
                          <input
                            type="checkbox"
                            checked={columnVisibility[columnId]}
                            disabled={isOnlyVisibleColumn}
                            onChange={() => handleColumnVisibilityChange(columnId)}
                          />
                          <span>{column.label}</span>
                        </label>
                        <div className="flex items-center gap-2">
                          <button
                            type="button"
                            className="rounded border border-gray-300 px-2 py-1 text-xs text-gray-700 disabled:cursor-not-allowed disabled:opacity-40"
                            disabled={index === 0}
                            onClick={() => moveColumn(columnId, 'left')}
                          >
                            Left
                          </button>
                          <button
                            type="button"
                            className="rounded border border-gray-300 px-2 py-1 text-xs text-gray-700 disabled:cursor-not-allowed disabled:opacity-40"
                            disabled={index === columnOrder.length - 1}
                            onClick={() => moveColumn(columnId, 'right')}
                          >
                            Right
                          </button>
                        </div>
                      </div>
                    );
                  })}
                </div>
                <div className="mt-4 flex items-center justify-between text-xs text-gray-500">
                  <span>Drag headers to reorder. Drag the header edge to resize. Use header filters and sort for visible columns.</span>
                  <button
                    type="button"
                    className="font-medium text-blue-600 hover:underline"
                    onClick={resetColumns}
                  >
                    Reset
                  </button>
                </div>
              </div>
            )}
          </div>
          <button
            type="button"
            onClick={() => onSelectVariants(visibleVariantIds)}
            disabled={visibleVariantIds.length === 0}
            className={`${TOOLBAR_BUTTON_CLASSES} font-medium text-slate-600 disabled:cursor-not-allowed disabled:opacity-50`}
          >
            Select all shown
          </button>
        </div>
        <div className="text-sm text-slate-500">
          {visibleColumns.length} of {RESULTS_TABLE_COLUMNS.length} columns shown
        </div>
      </div>
      {hasActiveTableFilters && (
        <div className="mb-3 text-sm text-gray-500">
          Table filters are active.
        </div>
      )}
      <div className="max-h-[75vh] overflow-auto rounded-xl border border-slate-200 bg-white">
        <table className="w-full border-collapse table-fixed" style={{ minWidth: Math.max(totalTableWidth, minimumTableWidth) }}>
          <colgroup>
            <col style={{ width: SELECTION_COLUMN_WIDTH }} />
            {visibleColumns.map((column) => {
              const isExpandingTxColumn = column.id === 'txImplications' && !hasExplicitTxColumnWidth;

              return (
                <col
                  key={column.id}
                  style={isExpandingTxColumn ? undefined : { width: columnWidths[column.id] }}
                />
              );
            })}
          </colgroup>
          <thead>
            <tr className="bg-slate-50">
              <th className="sticky top-0 z-10 border-r border-slate-100 bg-slate-50 p-0 shadow-[inset_0_-1px_0_0_rgb(226_232_240)]">
                <span className="sr-only">Select variant</span>
              </th>
              {visibleColumns.map((column) => (
                <th key={column.id} className="sticky top-0 z-10 border-r border-slate-100 bg-slate-50 p-0 text-left text-slate-500 shadow-[inset_0_-1px_0_0_rgb(226_232_240)] last:border-r-0">
                  <div
                    className="px-3 pb-2.5 pt-2.5 pr-5"
                  >
                    <div className="mb-1.5 flex items-start justify-between gap-2">
                      <div
                        className="cursor-grab break-words text-xs font-bold uppercase tracking-wide whitespace-normal"
                        draggable
                        onDragStart={() => setDraggedColumnId(column.id)}
                        onDragEnd={() => setDraggedColumnId(null)}
                        onDragOver={(event) => event.preventDefault()}
                        onDrop={() => handleDrop(column.id)}
                      >
                        {column.label}
                      </div>
                      {/* Sort arrows: dark for the sorted column (up = A-Z, down = Z-A), pale for the others */}
                      {(() => {
                        const activeSort = sortState ?? DEFAULT_SORT_STATE;
                        const isSorted = activeSort.columnId === column.id;
                        const SortIcon = !isSorted ? ArrowUpDown : activeSort.direction === 'asc' ? ArrowUp : ArrowDown;
                        const sortLabel = !isSorted ? `Sort by ${column.label}` : activeSort.direction === 'asc' ? `Sorted A-Z by ${column.label}` : `Sorted Z-A by ${column.label}`;
                        return (
                          <button
                            type="button"
                            aria-label={sortLabel}
                            title={sortLabel}
                            className={`shrink-0 rounded p-0.5 hover:bg-slate-200 ${isSorted ? 'text-slate-900' : 'text-slate-300 hover:text-slate-600'}`}
                            onClick={() => toggleSort(column.id)}
                          >
                            <SortIcon size={14} strokeWidth={2.5} aria-hidden="true" />
                          </button>
                        );
                      })()}
                    </div>
                    <input
                      type="text"
                      value={columnFilters[column.id]}
                      onChange={(event) => updateColumnFilter(column.id, event.target.value)}
                      placeholder={`Filter ${column.label.toLowerCase()}`}
                      className="w-full rounded-md border border-slate-200 bg-white px-2 py-1 text-xs font-normal text-slate-700 placeholder:text-slate-400"
                    />
                  </div>
                  <div
                    className="absolute right-0 top-0 h-full w-2 cursor-col-resize hover:bg-blue-200"
                    onMouseDown={(event) => startResize(event, column.id)}
                  />
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {processedGroups.map(([range, variantList]) => (
              <VariantGroupRow
                key={range}
                columns={visibleColumns}
                range={range}
                variants={variantList}
                renderDxImplications={renderDxImplications}
                renderTxImplications={renderTxImplications}
                renderMolecularConsequences={renderMolecularConsequences}
                oncogenicityResults={oncogenicityResults}
                getOncogenicityKey={getOncogenicityKey}
                onComputeOncogenicity={computeOncogenicityPrediction}
                onOpenOncogenicityDetails={openOncogenicityDetails}
                selectedVariantIds={selectedVariantIds}
                onToggleVariantSelection={onToggleVariantSelection}
                highlightedVariantId={highlightedVariantId}
              />
            ))}
            {processedGroups.length === 0 && (
              <tr>
                <td colSpan={visibleColumns.length + 1} className="px-6 py-12">
                  {/* Different messages before any search, after a search with no variants, and when table filters hide every row */}
                  <div className="flex flex-col items-center gap-2.5 text-center">
                    <span className={`flex h-11 w-11 items-center justify-center rounded-xl ${hasSearched && !hasActiveTableFilters ? 'bg-slate-100 text-slate-600' : 'bg-emerald-100 text-emerald-700'}`} aria-hidden="true">
                      {hasActiveTableFilters ? <Filter size={22} /> : hasSearched ? <Dna size={22} /> : <Search size={22} />}
                    </span>
                    <span className="text-base font-semibold text-slate-900">
                      {hasActiveTableFilters ? 'No rows match the table filters' : hasSearched ? 'No results found' : 'No search yet'}
                    </span>
                    <span className="max-w-md text-sm text-slate-500">
                      {hasActiveTableFilters
                        ? 'The filters under the column headers hide every row.'
                        : hasSearched
                          ? 'This patient has no variants in the searched genes or ranges.'
                          : 'Enter genes or genomic ranges above, or pick a preset, then click Search to see this patient\'s variants.'}
                    </span>
                    {hasActiveTableFilters && (
                      <button type="button" onClick={clearColumnFilters} className={TOOLBAR_BUTTON_CLASSES}>
                        Clear table filters
                      </button>
                    )}
                  </div>
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
      {selectedOncogenicityResult && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center bg-black/45 px-4 py-8"
          onClick={() => setSelectedOncogenicityKey(null)}
        >
          <div
            role="dialog"
            aria-modal="true"
            aria-labelledby="oncogenicity-modal-title"
            className="max-h-[85vh] w-full max-w-4xl overflow-y-auto rounded-2xl bg-white p-6 shadow-2xl"
            onClick={(event) => event.stopPropagation()}
          >
            <div className="mb-6 flex items-start justify-between gap-4">
              <div>
                <h2 id="oncogenicity-modal-title" className="text-2xl font-bold text-gray-900">
                  Oncogenicity Prediction
                </h2>
                <p className="mt-2 text-sm text-gray-600">
                  Review the overall prediction, evidence lines, and supporting caveats for this variant.
                </p>
                <p className="mt-2 text-sm text-gray-600">
                  Criteria and scoring rules follow the published oncogenicity recommendations described in{' '}
                  <a
                    href="https://pubmed.ncbi.nlm.nih.gov/35101336/"
                    target="_blank"
                    rel="noreferrer"
                    className="font-medium text-blue-700 underline decoration-blue-300 underline-offset-2 hover:text-blue-800"
                  >
                    Horak et al. 2022
                  </a>
                  , implemented using{' '}
                  <a
                    href="https://github.com/rhdolin/oncogenicity-predictor"
                    target="_blank"
                    rel="noreferrer"
                    className="font-medium text-blue-700 underline decoration-blue-300 underline-offset-2 hover:text-blue-800"
                  >
                    OncogenicityPredictor v1
                  </a>
                  .
                </p>
              </div>
              <button
                type="button"
                onClick={() => setSelectedOncogenicityKey(null)}
                className="rounded-md border border-gray-300 px-3 py-1.5 text-sm font-medium text-gray-700 transition-colors hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
              >
                Close
              </button>
            </div>

            <div className="mb-6 flex flex-wrap items-start gap-6 rounded-xl border border-gray-200 bg-gray-50 p-4">
              <div className="rounded-2xl bg-white p-3 shadow-sm ring-1 ring-slate-200">
                <Image
                  src={getGaugeImagePath(selectedOncogenicityResult.gauge)}
                  alt="Oncogenicity prediction gauge"
                  width={216}
                  height={122}
                  className="h-auto w-52"
                />
              </div>
              <dl className="grid min-w-0 flex-1 gap-3 sm:grid-cols-2">
                <div>
                  <dt className="text-xs font-semibold uppercase tracking-wide text-gray-500">Overall score</dt>
                  <dd className="mt-1 text-xl font-bold text-gray-900">
                    {selectedOncogenicityResult.score ?? 'Unavailable'}
                  </dd>
                </div>
                <div>
                  <dt className="text-xs font-semibold uppercase tracking-wide text-gray-500">Overall prediction</dt>
                  <dd className="mt-1 text-lg font-semibold text-gray-900">
                    {selectedOncogenicityResult.interpretation ?? 'Unavailable'}
                  </dd>
                </div>
                <div className="sm:col-span-2">
                  <dt className="text-xs font-semibold uppercase tracking-wide text-gray-500">Submitted Variant</dt>
                  <dd className="mt-1 break-all font-mono text-sm text-gray-700">
                    {selectedOncogenicityResult.variant}
                  </dd>
                </div>
              </dl>
            </div>

            <section className="mb-6">
              <h3 className="text-lg font-semibold text-gray-900">Evidence contributing to prediction</h3>
              <div className="mt-3 space-y-3">
                {selectedOncogenicityResult.observation?.component?.length ? selectedOncogenicityResult.observation.component.map((component, index) => (
                  <div key={`${component.code?.text ?? 'component'}-${index}`} className="rounded-xl border border-gray-200 bg-white p-4 shadow-sm">
                    <dl className="grid gap-3 sm:grid-cols-2">
                      <div>
                        <dt className="text-xs font-semibold uppercase tracking-wide text-gray-500">Evidence line</dt>
                        <dd className="mt-1 text-sm text-gray-900">{component.code?.text ?? 'Unavailable'}</dd>
                      </div>
                      <div>
                        <dt className="text-xs font-semibold uppercase tracking-wide text-gray-500">Score contribution</dt>
                        <dd className="mt-1 text-sm font-semibold text-gray-900">{component.valueInteger ?? 'Unavailable'}</dd>
                      </div>
                      <div>
                        <dt className="text-xs font-semibold uppercase tracking-wide text-gray-500">Criterion satisfied</dt>
                        <dd className="mt-1 text-sm text-gray-900">{getConceptDisplayText(component.interpretation) ?? 'Unavailable'}</dd>
                      </div>
                      <div>
                        <dt className="text-xs font-semibold uppercase tracking-wide text-gray-500">Details</dt>
                        <dd className="mt-1 text-sm text-gray-900 whitespace-normal break-words">{getConceptNarrativeText(component.interpretation) ?? 'Unavailable'}</dd>
                      </div>
                    </dl>
                  </div>
                )) : (
                  <div className="rounded-xl border border-dashed border-gray-300 bg-gray-50 p-4 text-sm text-gray-500">
                    No evidence lines were returned for this prediction.
                  </div>
                )}
              </div>
            </section>

            {getObservationCaveat(selectedOncogenicityResult.observation) && (
              <section className="mb-6 rounded-xl border border-amber-200 bg-amber-50 p-4">
                <h3 className="text-sm font-semibold uppercase tracking-wide text-amber-900">Caveat</h3>
                <p className="mt-2 whitespace-pre-wrap break-words text-sm text-amber-900">
                  {getObservationCaveat(selectedOncogenicityResult.observation)}
                </p>
              </section>
            )}

            <div className="flex flex-wrap items-center gap-3">
              <button
                type="button"
                onClick={loadExtendedEvidence}
                disabled={evidenceStatus === 'loading'}
                className={evidenceStatus === 'ready'
                  ? 'inline-flex items-center gap-2 rounded-md border border-gray-300 bg-gray-100 px-4 py-2.5 text-sm font-semibold text-gray-500 transition-colors hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2'
                  : 'inline-flex items-center gap-2 rounded-md border border-blue-700 bg-blue-600 px-4 py-2.5 text-sm font-semibold text-white shadow-[0_3px_0_0_rgb(29_78_216)] transition-[transform,box-shadow,background-color] hover:bg-blue-700 hover:shadow-[0_2px_0_0_rgb(30_64_175)] active:translate-y-px active:shadow-[0_1px_0_0_rgb(30_64_175)] focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2 disabled:cursor-not-allowed disabled:opacity-60'}
              >
                {evidenceStatus === 'loading' && (
                  <span className="h-4 w-4 flex-shrink-0 animate-spin rounded-full border-2 border-blue-200 border-t-white" />
                )}
                {evidenceButtonLabel}
                {evidenceStatus === 'ready' && <ArrowDownToLine className="h-4 w-4" aria-hidden="true" />}
              </button>
              <button
                type="button"
                onClick={() => setShowClinVarWorkInProgress((currentValue) => !currentValue)}
                className="rounded-md border border-gray-300 px-4 py-2.5 text-sm font-semibold text-gray-700 transition-colors hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-blue-300 focus:ring-offset-2"
              >
                Submit to ClinVar
              </button>
            </div>

            {showClinVarWorkInProgress && (
              <div className="mt-3 rounded-lg border border-gray-200 bg-gray-50 px-4 py-3 text-sm text-gray-700">
                Work in progress...
              </div>
            )}

            {selectedOncogenicityResult.evidenceStatus === 'error' && selectedOncogenicityResult.evidenceError && (
              <div className="mt-4 rounded-lg border border-red-200 bg-red-50 px-4 py-3 text-sm text-red-700">
                {selectedOncogenicityResult.evidenceError}
              </div>
            )}

            {selectedOncogenicityResult.evidenceStatus === 'ready' && selectedOncogenicityResult.evidenceJson !== undefined && (
              <section ref={extendedEvidenceRef} className="mt-6 scroll-mt-4">
                <h3 className="text-lg font-semibold text-gray-900">Extended evidence details</h3>
                <pre className="mt-3 overflow-x-auto rounded-xl border border-gray-200 bg-gray-950 p-4 text-xs leading-6 text-gray-100 whitespace-pre-wrap break-all">
                  {JSON.stringify(selectedOncogenicityResult.evidenceJson, null, 2)}
                </pre>
              </section>
            )}
          </div>
        </div>
      )}
    </div>
  );
}
