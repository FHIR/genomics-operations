import { useState, useMemo, useEffect } from 'react';
import { Sidebar } from '@/components/sidebar';
import { Variant } from '@/types/variants';

export type FilterSource = 'user' | 'radio' | 'system';

export interface ActiveFilter<T> {
  value: T;
  source: FilterSource;
}

export interface FilterCriteria {
  selectedImpacts: string[];
  selectedConsequences: string[];
  selectedDxSignificances: string[];
  selectedDxStars: number[];
  selectedTxEvidence: ActiveFilter<string>[];
  selectedTxMedications: string[];
  selectedTxImplications: string[];
  selectedTxPhenotypes: string[];
  actionableFilter?: string;
  evidenceLevelSource?: FilterSource;
}

interface FilterSidebarProps {
  isOpen: boolean;
  onClose: () => void;
  onFilterChange: (filters: FilterCriteria) => void;
  results: Variant[];
  currentFilters?: FilterCriteria;
}

interface FilterSectionProps {
  title: string;
  isOpen: boolean;
  onToggle: () => void;
  children: React.ReactNode;
}

function FilterSection({ title, isOpen, onToggle, children }: FilterSectionProps) {
  return (
    <div className="border rounded-lg mb-2">
      <button
        onClick={onToggle}
        className="w-full px-4 py-2 flex items-center justify-between bg-gray-50 hover:bg-gray-100 rounded-t-lg"
      >
        <span className="font-medium">{title}</span>
        <svg
          xmlns="http://www.w3.org/2000/svg"
          className={`h-5 w-5 transform transition-transform ${isOpen ? 'rotate-180' : ''}`}
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>
      {isOpen && (
        <div className="p-4 border-t">
          {children}
        </div>
      )}
    </div>
  );
}

export default function FilterSidebar({ isOpen, onClose, onFilterChange, results, currentFilters }: FilterSidebarProps) {
  const [filters, setFilters] = useState<FilterCriteria>(currentFilters || {
    selectedImpacts: [],
    selectedConsequences: [],
    selectedDxSignificances: [],
    selectedDxStars: [],
    selectedTxEvidence: [],
    selectedTxMedications: [],
    selectedTxImplications: [],
    selectedTxPhenotypes: []
  });

  // Sync filters with currentFilters prop
  useEffect(() => {
    if (currentFilters) {
      setFilters(currentFilters);
    }
  }, [currentFilters]);

  const [openSections, setOpenSections] = useState({
    molecular: false,
    therapeutic: false,
    dx: false
  });

  const toggleSection = (section: keyof typeof openSections) => {
    setOpenSections(prev => ({ ...prev, [section]: !prev[section] }));
  };

  // Extract unique impact levels dynamically from results
  const impactLevels = useMemo(() => {
    const impactSet = new Set<string>();

    results.forEach(variant => {
      if (variant.molecularConsequences) {
        variant.molecularConsequences.forEach(mc => {
          if (mc.impact) {
            impactSet.add(mc.impact);
          }
        });
      }
    });

    // Convert to array and sort by severity (HIGH -> MODERATE -> LOW -> MODIFIER -> others)
    const severityOrder = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER'];
    return Array.from(impactSet).sort((a, b) => {
      const aIndex = severityOrder.indexOf(a);
      const bIndex = severityOrder.indexOf(b);

      // If both are in the severity order, sort by that order
      if (aIndex !== -1 && bIndex !== -1) {
        return aIndex - bIndex;
      }
      // If only one is in the severity order, prioritize it
      if (aIndex !== -1) return -1;
      if (bIndex !== -1) return 1;
      // If neither is in the severity order, sort alphabetically
      return a.localeCompare(b);
    });
  }, [results]);

  // Extract unique molecular consequences dynamically from results
  const molecularImplications = useMemo(() => {
    const consequenceSet = new Set<string>();

    results.forEach(variant => {
      if (variant.molecularConsequences) {
        variant.molecularConsequences.forEach(mc => {
          if (mc.featureConsequence) {
            consequenceSet.add(mc.featureConsequence);
          }
        });
      }
    });

    // Convert to array and sort alphabetically
    return Array.from(consequenceSet).sort();
  }, [results]);

  // Extract unique DX significances dynamically from results
  const dxSignificances = useMemo(() => {
    const significanceSet = new Set<string>();

    results.forEach(variant => {
      if (variant.dxImplications) {
        variant.dxImplications.forEach(dx => {
          if (dx.clinicalSignificance) {
            significanceSet.add(dx.clinicalSignificance);
          }
        });
      }
    });

    // Convert to array and sort alphabetically
    return Array.from(significanceSet).sort();
  }, [results]);

  // Extract unique DX star ratings dynamically from results
  const dxStars = useMemo(() => {
    const starsSet = new Set<number>();

    results.forEach(variant => {
      if (variant.dxImplications) {
        variant.dxImplications.forEach(dx => {
          if (dx.evidenceLevel) {
            // Calculate stars for this implication
            const evidenceLevel = dx.evidenceLevel || '';
            let stars;
            switch (evidenceLevel.toLowerCase()) {
              case 'practice guideline':
                stars = 4;
                break;
              case 'reviewed by expert panel':
                stars = 3;
                break;
              case 'criteria provided, single submitter':
                stars = 1;
                break;
              case 'no assertion criteria provided':
              case 'no classification provided':
              default:
                stars = 0;
                break;
            }
            starsSet.add(stars);
          }
        });
      }
    });

    // Convert to array and sort numerically (0, 1, 3, 4)
    return Array.from(starsSet).sort((a, b) => a - b);
  }, [results]);

  // Extract unique therapeutic evidence levels dynamically from results
  const txEvidenceLevels = useMemo(() => {
    const evidenceSet = new Set<string>();

    results.forEach(variant => {
      if (variant.txImplications) {
        variant.txImplications.forEach(tx => {
          if (tx.evidenceLevel) {
            evidenceSet.add(tx.evidenceLevel);
          }
        });
      }
    });

    return Array.from(evidenceSet).sort();
  }, [results]);

  // Extract unique therapeutic medications dynamically from results
  const txMedications = useMemo(() => {
    const medicationSet = new Set<string>();

    results.forEach(variant => {
      if (variant.txImplications) {
        variant.txImplications.forEach(tx => {
          if (tx.medication && tx.medication.trim() !== '') {
            medicationSet.add(tx.medication);
          }
        });
      }
    });

    return Array.from(medicationSet).sort();
  }, [results]);

  // Extract unique therapeutic implications dynamically from results
  const txImplications = useMemo(() => {
    const implicationSet = new Set<string>();

    results.forEach(variant => {
      if (variant.txImplications) {
        variant.txImplications.forEach(tx => {
          if (tx.therapeuticImplication && tx.therapeuticImplication.trim() !== '') {
            implicationSet.add(tx.therapeuticImplication);
          }
        });
      }
    });

    return Array.from(implicationSet).sort();
  }, [results]);

  // Extract unique therapeutic phenotypes dynamically from results
  const txPhenotypes = useMemo(() => {
    const phenotypeSet = new Set<string>();

    results.forEach(variant => {
      if (variant.txImplications) {
        variant.txImplications.forEach(tx => {
          if (tx.phenotypicContext) {
            phenotypeSet.add(tx.phenotypicContext);
          }
        });
      }
    });

    return Array.from(phenotypeSet).sort();
  }, [results]);

  const handleSelectedImpactChange = (impact: string) => {
    const newSelectedImpacts = filters.selectedImpacts.includes(impact)
      ? filters.selectedImpacts.filter(i => i !== impact)
      : [...filters.selectedImpacts, impact];

    const newFilters = { ...filters, selectedImpacts: newSelectedImpacts };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleSelectedConsequenceChange = (consequence: string) => {
    const newSelectedConsequences = filters.selectedConsequences.includes(consequence)
      ? filters.selectedConsequences.filter(c => c !== consequence)
      : [...filters.selectedConsequences, consequence];

    const newFilters = { ...filters, selectedConsequences: newSelectedConsequences };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleSelectedDxSignificanceChange = (significance: string) => {
    const newSelectedDxSignificances = filters.selectedDxSignificances.includes(significance)
      ? filters.selectedDxSignificances.filter(s => s !== significance)
      : [...filters.selectedDxSignificances, significance];

    const newFilters = { ...filters, selectedDxSignificances: newSelectedDxSignificances };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleSelectedDxStarChange = (stars: number) => {
    const newSelectedDxStars = filters.selectedDxStars.includes(stars)
      ? filters.selectedDxStars.filter(s => s !== stars)
      : [...filters.selectedDxStars, stars];

    const newFilters = { ...filters, selectedDxStars: newSelectedDxStars };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleSelectedTxEvidenceChange = (evidence: string) => {
    const existingFilterIndex = filters.selectedTxEvidence.findIndex(filter => filter.value === evidence);
    let newSelectedTxEvidence: ActiveFilter<string>[];

    if (existingFilterIndex !== -1) {
      // Remove the filter if it exists and is user-added
      const existingFilter = filters.selectedTxEvidence[existingFilterIndex];
      if (existingFilter.source === 'user') {
        newSelectedTxEvidence = filters.selectedTxEvidence.filter(f => f.value !== evidence);
      } else {
        // Can't remove system/radio filters
        return;
      }
    } else {
      // Add as user filter
      newSelectedTxEvidence = [...filters.selectedTxEvidence, { value: evidence, source: 'user' }];
    }

    const newFilters = { ...filters, selectedTxEvidence: newSelectedTxEvidence };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleSelectedTxMedicationChange = (medication: string) => {
    const newSelectedTxMedications = filters.selectedTxMedications.includes(medication)
      ? filters.selectedTxMedications.filter(m => m !== medication)
      : [...filters.selectedTxMedications, medication];

    const newFilters = { ...filters, selectedTxMedications: newSelectedTxMedications };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleSelectedTxImplicationChange = (implication: string) => {
    const newSelectedTxImplications = filters.selectedTxImplications.includes(implication)
      ? filters.selectedTxImplications.filter(i => i !== implication)
      : [...filters.selectedTxImplications, implication];

    const newFilters = { ...filters, selectedTxImplications: newSelectedTxImplications };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleSelectedTxPhenotypeChange = (phenotype: string) => {
    const newSelectedTxPhenotypes = filters.selectedTxPhenotypes.includes(phenotype)
      ? filters.selectedTxPhenotypes.filter(p => p !== phenotype)
      : [...filters.selectedTxPhenotypes, phenotype];

    const newFilters = { ...filters, selectedTxPhenotypes: newSelectedTxPhenotypes };
    setFilters(newFilters);
    onFilterChange(newFilters);
  };

  const handleReset = () => {
    // Preserve radio button filters when resetting
    const radioFilters = filters.selectedTxEvidence.filter(f => f.source === 'radio');

    const resetFilters: FilterCriteria = {
      selectedImpacts: [],
      selectedConsequences: [],
      selectedDxSignificances: [],
      selectedDxStars: [],
      selectedTxEvidence: radioFilters, // Keep radio button filters
      selectedTxMedications: [],
      selectedTxImplications: [],
      selectedTxPhenotypes: [],
      // Preserve the actionable filter info
      actionableFilter: filters.actionableFilter,
      evidenceLevelSource: filters.evidenceLevelSource
    };
    setFilters(resetFilters);
    onFilterChange(resetFilters);
  };

  return (
    <Sidebar
      isOpen={isOpen}
      onClose={onClose}
      title="Filters"
      width="sm"
    >
      {/* Scrollable content area */}
      <div className="flex-1 overflow-y-auto">
        <div className="p-4">
          {/* Actionable Filter Status Section */}
          {filters.actionableFilter && (
            <div className="mb-4 p-3 bg-blue-50 border border-blue-200 rounded-lg">
              <h4 className="font-medium text-blue-900 mb-2">Active Radio Button Filter</h4>
              <div className="text-sm text-blue-800">
                <div className="mb-1">
                  <span className="font-medium">Filter Type:</span> {filters.actionableFilter}
                </div>
                {filters.selectedTxEvidence.filter(f => f.source === 'radio').length > 0 && (
                  <div className="mb-1">
                    <span className="font-medium">Applied Evidence Levels:</span>{' '}
                    {filters.selectedTxEvidence
                      .filter(f => f.source === 'radio')
                      .map(f => f.value)
                      .join(', ')}
                  </div>
                )}
                <div className="text-xs text-blue-600">
                  {filters.actionableFilter === "Possibly actionable"
                    ? "This filter selects all non-A evidence levels (B, C, D, E, etc.) and excludes A-level evidence."
                    : "This filter requires A-level evidence and applies phenotype matching based on your selected cancer type."
                  }
                </div>
              </div>
            </div>
          )}

          <FilterSection
            title="Molecular Consequences"
            isOpen={openSections.molecular}
            onToggle={() => toggleSection('molecular')}
          >
            <div className="space-y-4">
              <div>
                <h4 className="font-medium mb-2 mt-4">Impact Level</h4>
                <div className="space-y-2">
                  {impactLevels.map(impact => (
                    <label key={impact} className="flex items-center">
                      <input
                        type="checkbox"
                        checked={filters.selectedImpacts.includes(impact)}
                        onChange={() => handleSelectedImpactChange(impact)}
                        className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                      />
                       <span className="text-sm" style={{ marginLeft: '4px' }}>{impact}</span>
                    </label>
                  ))}
                </div>
              </div>

              <div>
                <h4 className="font-medium mb-2 mt-4">Consequence Type</h4>
                <div className="space-y-2 max-h-60 overflow-y-auto">
                  {molecularImplications.map(consequence => (
                    <label key={consequence} className="flex items-center">
                      <input
                        type="checkbox"
                        checked={filters.selectedConsequences.includes(consequence)}
                        onChange={() => handleSelectedConsequenceChange(consequence)}
                        className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                      />
                      <span className="text-sm" style={{ marginLeft: '4px' }}>{consequence}</span>
                    </label>
                  ))}
                </div>
              </div>
            </div>
          </FilterSection>

          <FilterSection
            title="Therapeutic Implications"
            isOpen={openSections.therapeutic}
            onToggle={() => toggleSection('therapeutic')}
          >
            <div className="space-y-4">
              <div>
                <h4 className="font-medium mb-2 mt-4">Evidence Level</h4>
                <div className="space-y-2 max-h-40 overflow-y-auto">
                  {txEvidenceLevels.map(evidence => {
                    const existingFilter = filters.selectedTxEvidence.find(f => f.value === evidence);
                    const isChecked = !!existingFilter;
                    const isDisabled = existingFilter?.source === 'radio' || existingFilter?.source === 'system';

                    return (
                      <label key={evidence} className={`flex items-center ${isDisabled ? 'opacity-75' : ''}`}>
                        <input
                          type="checkbox"
                          checked={isChecked}
                          disabled={isDisabled}
                          onChange={() => handleSelectedTxEvidenceChange(evidence)}
                          className={`rounded border-gray-300 focus:ring-blue-500 ${
                            isDisabled ? 'text-gray-400' : 'text-blue-600'
                          }`}
                        />
                        <span className="text-sm" style={{ marginLeft: '4px' }}>
                          {evidence}
                          {existingFilter?.source === 'radio' && (
                            <span className="ml-2 px-2 py-1 bg-blue-100 text-blue-800 text-xs rounded">
                              Actionable Filter
                            </span>
                          )}
                          {existingFilter?.source === 'system' && (
                            <span className="ml-2 px-2 py-1 bg-gray-100 text-gray-800 text-xs rounded">
                              System
                            </span>
                          )}
                        </span>
                      </label>
                    );
                  })}
                </div>
              </div>

              {txMedications.length > 0 && (
                <div>
                  <h4 className="font-medium mb-2 mt-4">Medication</h4>
                  <div className="space-y-2 max-h-40 overflow-y-auto">
                    {txMedications.map(medication => (
                      <label key={medication} className="flex items-center">
                        <input
                          type="checkbox"
                          checked={filters.selectedTxMedications.includes(medication)}
                          onChange={() => handleSelectedTxMedicationChange(medication)}
                          className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                        />
                        <span className="text-sm" style={{ marginLeft: '4px' }}>{medication}</span>
                      </label>
                    ))}
                  </div>
                </div>
              )}

              {txImplications.length > 0 && (
                <div>
                  <h4 className="font-medium mb-2 mt-4">Therapeutic Implication</h4>
                  <div className="space-y-2 max-h-40 overflow-y-auto">
                    {txImplications.map(implication => (
                      <label key={implication} className="flex items-center">
                        <input
                          type="checkbox"
                          checked={filters.selectedTxImplications.includes(implication)}
                          onChange={() => handleSelectedTxImplicationChange(implication)}
                          className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                        />
                        <span className="text-sm" style={{ marginLeft: '4px' }}>{implication}</span>
                      </label>
                    ))}
                  </div>
                </div>
              )}

              <div>
                <h4 className="font-medium mb-2 mt-4">Phenotypic Context</h4>
                <div className="space-y-2 max-h-40 overflow-y-auto">
                  {txPhenotypes.map(phenotype => (
                    <label key={phenotype} className="flex items-center">
                      <input
                        type="checkbox"
                        checked={filters.selectedTxPhenotypes.includes(phenotype)}
                        onChange={() => handleSelectedTxPhenotypeChange(phenotype)}
                        className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                      />
                      <span className="text-sm" style={{ marginLeft: '4px' }}>{phenotype}</span>
                    </label>
                  ))}
                </div>
              </div>
            </div>
          </FilterSection>

          <FilterSection
            title="Diagnostic Implications"
            isOpen={openSections.dx}
            onToggle={() => toggleSection('dx')}
          >
            <div className="space-y-4">
              <div>
                <h4 className="font-medium mb-2 mt-4">Clinical Significance</h4>
                <div className="space-y-2 max-h-40 overflow-y-auto">
                  {dxSignificances.map(significance => (
                    <label key={significance} className="flex items-center">
                      <input
                        type="checkbox"
                        checked={filters.selectedDxSignificances.includes(significance)}
                        onChange={() => handleSelectedDxSignificanceChange(significance)}
                        className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                      />
                      <span className="text-sm" style={{ marginLeft: '4px' }}>{significance}</span>
                    </label>
                  ))}
                </div>
              </div>

              <div>
                <h4 className="font-medium mb-2 mt-4">ClinVar Stars</h4>
                <div className="space-y-2">
                  {dxStars.map(stars => (
                    <label key={stars} className="flex items-center">
                      <input
                        type="checkbox"
                        checked={filters.selectedDxStars.includes(stars)}
                        onChange={() => handleSelectedDxStarChange(stars)}
                        className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                      />
                      <span className="text-sm" style={{ marginLeft: '4px' }}>{stars} {stars === 1 ? 'star' : 'stars'}</span>
                    </label>
                  ))}
                </div>
              </div>
            </div>
          </FilterSection>
        </div>
      </div>

      {/* Fixed footer */}
      <div className="p-4 border-t bg-white flex-shrink-0">
        <div className="flex justify-between items-center">
          <button
            onClick={handleReset}
            className="px-4 py-2 bg-gray-100 text-gray-700 rounded hover:bg-gray-200 transition-colors"
          >
            Reset Filters
          </button>
          <button
            onClick={onClose}
            className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 transition-colors"
          >
            Close
          </button>
        </div>
        {filters.actionableFilter && (
          <div className="mt-2 text-xs text-gray-600 text-center">
            To clear radio button filters, change the selection above the search box
          </div>
        )}
      </div>
    </Sidebar>
  );
}