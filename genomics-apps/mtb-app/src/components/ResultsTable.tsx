"use client";

import { MouseEvent as ReactMouseEvent, useEffect, useMemo, useState } from 'react';
import { Variant } from '@/types/variants';
import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import DxImplicationCell from './DxImplicationCell';
import MolecularConsequenceCell from './MolecularConsequenceCell';
import TxImplicationCell from './TxImplicationCell';
import VariantGroupRow from './VariantGroupRow';
import {
  ColumnId,
  DEFAULT_COLUMN_ORDER,
  RESULTS_TABLE_COLUMNS,
} from './resultsTableColumns';
import _ from 'lodash';

type SortDirection = 'asc' | 'desc';

interface SortState {
  columnId: ColumnId;
  direction: SortDirection;
}

const COLUMN_MAP = Object.fromEntries(
  RESULTS_TABLE_COLUMNS.map((column) => [column.id, column])
) as Record<ColumnId, (typeof RESULTS_TABLE_COLUMNS)[number]>;

const DEFAULT_COLUMN_VISIBILITY = Object.fromEntries(
  RESULTS_TABLE_COLUMNS.map((column) => [column.id, column.id !== 'oncogenicityPrediction'])
) as Record<ColumnId, boolean>;

const DEFAULT_COLUMN_WIDTHS = Object.fromEntries(
  RESULTS_TABLE_COLUMNS.map((column) => [column.id, column.defaultWidth])
) as Record<ColumnId, number>;

const TABLE_SETTINGS_STORAGE_KEY = 'mtb-results-table-settings';

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
  onToggleFilters: () => void;
}

function getColumnTextValue(variant: Variant, columnId: ColumnId, range: string) {
  switch (columnId) {
    case 'range':
      return range;
    case 'variant':
      return [variant.variant, variant.molecularConsequences?.[0]?.proteinChange]
        .filter(Boolean)
        .join(' ');
    case 'oncogenicityPrediction':
      return variant.oncogenicityPrediction ?? '';
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

  const comparison = normalizedLeft.localeCompare(normalizedRight, undefined, {
    numeric: true,
    sensitivity: 'base',
  });

  return direction === 'asc' ? comparison : -comparison;
}

export default function ResultsTable({ results, onToggleFilters }: ResultsTableProps) {
  const [columnOrder, setColumnOrder] = useState<ColumnId[]>(DEFAULT_COLUMN_ORDER);
  const [columnVisibility, setColumnVisibility] = useState<Record<ColumnId, boolean>>(
    DEFAULT_COLUMN_VISIBILITY
  );
  const [columnWidths, setColumnWidths] = useState<Record<ColumnId, number>>(
    DEFAULT_COLUMN_WIDTHS
  );
  const [draggedColumnId, setDraggedColumnId] = useState<ColumnId | null>(null);
  const [sortState, setSortState] = useState<SortState | null>(null);
  const [columnFilters, setColumnFilters] = useState<Record<ColumnId, string>>({
    range: '',
    variant: '',
    oncogenicityPrediction: '',
    molecularConsequences: '',
    dxImplications: '',
    txImplications: '',
  });

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
        columnWidths?: Partial<Record<ColumnId, number>>;
      };

      setColumnOrder(normalizeColumnOrder(parsedSettings.columnOrder));
      setColumnVisibility({
        ...DEFAULT_COLUMN_VISIBILITY,
        ...parsedSettings.columnVisibility,
      });
      setColumnWidths({
        ...DEFAULT_COLUMN_WIDTHS,
        ...Object.fromEntries(
          Object.entries(parsedSettings.columnWidths ?? {}).map(([columnId, width]) => {
            const typedColumnId = columnId as ColumnId;
            const minimumWidth = COLUMN_MAP[typedColumnId]?.minWidth ?? 140;

            return [typedColumnId, Math.max(Number(width) || minimumWidth, minimumWidth)];
          })
        ),
      });
    } catch {
      window.localStorage.removeItem(TABLE_SETTINGS_STORAGE_KEY);
    }
  }, []);

  useEffect(() => {
    if (typeof window === 'undefined') {
      return;
    }

    window.localStorage.setItem(
      TABLE_SETTINGS_STORAGE_KEY,
      JSON.stringify({
        columnOrder,
        columnVisibility,
        columnWidths,
      })
    );
  }, [columnOrder, columnVisibility, columnWidths]);

  // Function to render Dx Implications content
  const renderDxImplications = (implications?: DxImplication[]) => {
    return <DxImplicationCell implications={implications} />;
  };

  // Function to render Molecular Consequences content
  const renderMolecularConsequences = (consequences?: MolecularConsequence[]) => {
    return <MolecularConsequenceCell consequences={consequences} />;
  };

  // Function to render Tx Implications content
  const renderTxImplications = (implications?: ProcessedTxImplication[]) => {
    return <TxImplicationCell implications={implications} />;
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
    const activeVisibleFilters = visibleColumns.filter(
      (column) => columnFilters[column.id].trim() !== ''
    );

    const filteredGroups = groupedResults
      .map(([range, variants]) => {
        const filteredVariants = variants.filter((variant) => {
          return activeVisibleFilters.every((column) => {
            const filterValue = columnFilters[column.id].trim().toLocaleLowerCase();
            const cellValue = getColumnTextValue(variant, column.id, range).toLocaleLowerCase();
            return cellValue.includes(filterValue);
          });
        });

        return [range, filteredVariants] as const;
      })
      .filter(([, variants]) => variants.length > 0);

    if (!sortState || !columnVisibility[sortState.columnId]) {
      return filteredGroups;
    }

    const sortedGroups = filteredGroups
      .map(([range, variants]) => {
        const sortedVariants = [...variants].sort((leftVariant, rightVariant) => {
          const leftValue = getColumnTextValue(leftVariant, sortState.columnId, range);
          const rightValue = getColumnTextValue(rightVariant, sortState.columnId, range);
          return compareColumnValues(leftValue, rightValue, sortState.direction);
        });

        return [range, sortedVariants] as const;
      })
      .sort(([leftRange, leftVariants], [rightRange, rightVariants]) => {
        const leftValue = getColumnTextValue(leftVariants[0], sortState.columnId, leftRange);
        const rightValue = getColumnTextValue(rightVariants[0], sortState.columnId, rightRange);
        return compareColumnValues(leftValue, rightValue, sortState.direction);
      });

    return sortedGroups;
  }, [groupedResults, visibleColumns, columnFilters, sortState, columnVisibility]);

  const totalTableWidth = visibleColumns.reduce(
    (width, column) => width + (columnWidths[column.id] ?? column.defaultWidth),
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
    setSortState(null);
    setColumnFilters({
      range: '',
      variant: '',
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
    const minimumWidth = COLUMN_MAP[columnId].minWidth;

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

  return (
    <div className="bg-gray-100 p-6 rounded-lg shadow-sm w-full min-w-[1400px] -ml-64">
      <div className="mb-4 flex flex-wrap items-center justify-between gap-3 rounded-lg border border-gray-200 bg-white px-4 py-3">
        <div className="flex flex-wrap items-center gap-3">
          <button
            type="button"
            onClick={onToggleFilters}
            className="px-4 py-2 bg-white text-gray-700 rounded hover:bg-gray-50 transition-colors flex items-center gap-2 border border-gray-300"
          >
            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5 text-gray-500" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
            </svg>
            Filter Results
          </button>
          <details className="relative">
            <summary className="list-none px-4 py-2 bg-white text-gray-700 rounded hover:bg-gray-50 transition-colors flex items-center gap-2 border border-gray-300 cursor-pointer">
              Customize Table
            </summary>
            <div className="absolute left-0 z-10 mt-2 w-96 rounded-lg border border-gray-200 bg-white p-4 shadow-xl">
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
          </details>
        </div>
        <div className="text-sm text-gray-500">
          {visibleColumns.length} of {RESULTS_TABLE_COLUMNS.length} columns shown
        </div>
      </div>
      {hasActiveTableFilters && (
        <div className="mb-3 text-sm text-gray-500">
          Table filters are active.
        </div>
      )}
      <div className="overflow-x-auto">
        <table className="border-collapse table-fixed" style={{ minWidth: totalTableWidth }}>
          <colgroup>
            {visibleColumns.map((column) => (
              <col key={column.id} style={{ width: columnWidths[column.id] }} />
            ))}
          </colgroup>
          <thead>
            <tr className="bg-gray-200">
              {visibleColumns.map((column) => (
                <th key={column.id} className="relative border-r border-gray-300 p-0 text-left text-gray-700 last:border-r-0">
                  <div
                    className="px-3 pb-3 pt-2 pr-5"
                  >
                    <div className="mb-2 flex items-start justify-between gap-2">
                      <div
                        className="cursor-grab text-sm font-semibold"
                        draggable
                        onDragStart={() => setDraggedColumnId(column.id)}
                        onDragEnd={() => setDraggedColumnId(null)}
                        onDragOver={(event) => event.preventDefault()}
                        onDrop={() => handleDrop(column.id)}
                      >
                        {column.label}
                      </div>
                      <button
                        type="button"
                        className={`rounded border px-2 py-1 text-xs ${sortState?.columnId === column.id ? 'border-blue-500 bg-blue-50 text-blue-700' : 'border-gray-300 bg-white text-gray-700'}`}
                        onClick={() => toggleSort(column.id)}
                      >
                        {sortState?.columnId === column.id
                          ? sortState.direction === 'asc'
                            ? 'Sort A-Z'
                            : 'Sort Z-A'
                          : 'Sort'}
                      </button>
                    </div>
                    <input
                      type="text"
                      value={columnFilters[column.id]}
                      onChange={(event) => updateColumnFilter(column.id, event.target.value)}
                      placeholder={`Filter ${column.label.toLowerCase()}`}
                      className="w-full rounded border border-gray-300 px-2 py-1 text-xs font-normal text-gray-700 placeholder:text-gray-400"
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
              />
            ))}
            {processedGroups.length === 0 && (
              <tr>
                <td colSpan={visibleColumns.length} className="p-4 text-center text-gray-500">
                  {hasActiveTableFilters ? 'No rows match the current table filters' : 'No results found'}
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
}
