import { FilterCriteria } from './FilterSidebar';
import { Variant } from '@/types/variants';

export function applyFiltersToVariants(variants: Variant[], filters: FilterCriteria): Variant[] {
    return variants.filter(variant => {
        // Filter by molecular consequences
        const hasMatchingConsequence = filters.selectedConsequences.length === 0 ||
            variant.molecularConsequences?.some(mc =>
                filters.selectedConsequences.includes(mc.featureConsequence || ''));

        const hasMatchingImpact = filters.selectedImpacts.length === 0 ||
            variant.molecularConsequences?.some(mc =>
                filters.selectedImpacts.includes(mc.impact || ''));

        // Filter by therapeutic implications
        const hasMatchingTxImplication = (
            filters.selectedTxEvidence.length === 0 &&
            filters.selectedTxMedications.length === 0 &&
            filters.selectedTxImplications.length === 0 &&
            filters.selectedTxPhenotypes.length === 0
        ) || variant.txImplications?.some(tx => {
            const matchesEvidence = filters.selectedTxEvidence.length === 0 ||
                filters.selectedTxEvidence.some(filter =>
                    tx.evidenceLevel && filter.value === tx.evidenceLevel);

            const matchesMedication = filters.selectedTxMedications.length === 0 ||
                (tx.medication && tx.medication.trim() !== '' && filters.selectedTxMedications.includes(tx.medication));

            const matchesImplication = filters.selectedTxImplications.length === 0 ||
                (tx.therapeuticImplication && tx.therapeuticImplication.trim() !== '' && filters.selectedTxImplications.includes(tx.therapeuticImplication));

            const matchesPhenotype = filters.selectedTxPhenotypes.length === 0 ||
                (tx.phenotypicContext && filters.selectedTxPhenotypes.includes(tx.phenotypicContext));

            return matchesEvidence && matchesMedication && matchesImplication && matchesPhenotype;
        });

        // Filter by DX significance
        const hasMatchingDxSignificance = filters.selectedDxSignificances.length === 0 ||
            variant.dxImplications?.some(dx =>
                dx.clinicalSignificance && filters.selectedDxSignificances.includes(dx.clinicalSignificance));

        // Filter by DX stars
        const hasMatchingDxStars = filters.selectedDxStars.length === 0 ||
            variant.dxImplications?.some(dx => {
                if (!dx.evidenceLevel) return false;

                // Calculate stars for this implication (same logic as in FilterSidebar)
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

                return filters.selectedDxStars.includes(stars);
            });

        return hasMatchingConsequence &&
               hasMatchingImpact &&
               hasMatchingTxImplication &&
               hasMatchingDxSignificance &&
               hasMatchingDxStars;
    }).map(variant => {
        // Create a copy of the variant with filtered implications
        const filteredVariant = { ...variant };

        // Filter therapeutic implications if there are TX filters active
        if (filters.selectedTxEvidence.length > 0 ||
            filters.selectedTxMedications.length > 0 ||
            filters.selectedTxImplications.length > 0 ||
            filters.selectedTxPhenotypes.length > 0) {
            filteredVariant.txImplications = variant.txImplications?.filter(tx => {
                const matchesEvidence = filters.selectedTxEvidence.length === 0 ||
                    filters.selectedTxEvidence.some(filter =>
                        tx.evidenceLevel && filter.value === tx.evidenceLevel);

                const matchesMedication = filters.selectedTxMedications.length === 0 ||
                    (tx.medication && tx.medication.trim() !== '' && filters.selectedTxMedications.includes(tx.medication));

                const matchesImplication = filters.selectedTxImplications.length === 0 ||
                    (tx.therapeuticImplication && tx.therapeuticImplication.trim() !== '' && filters.selectedTxImplications.includes(tx.therapeuticImplication));

                const matchesPhenotype = filters.selectedTxPhenotypes.length === 0 ||
                    (tx.phenotypicContext && filters.selectedTxPhenotypes.includes(tx.phenotypicContext));

                return matchesEvidence && matchesMedication && matchesImplication && matchesPhenotype;
            });
        }

        // Filter diagnostic implications if there are DX filters active
        if (filters.selectedDxSignificances.length > 0 || filters.selectedDxStars.length > 0) {
            filteredVariant.dxImplications = variant.dxImplications?.filter(dx => {
                // Check significance filter
                const matchesSignificance = filters.selectedDxSignificances.length === 0 ||
                    (dx.clinicalSignificance && filters.selectedDxSignificances.includes(dx.clinicalSignificance));

                // Check stars filter
                const matchesStars = filters.selectedDxStars.length === 0 ||
                    (() => {
                        if (!dx.evidenceLevel) return false;

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

                        return filters.selectedDxStars.includes(stars);
                    })();

                return matchesSignificance && matchesStars;
            });
        }

        // Filter molecular consequences if there are MC filters active
        if (filters.selectedConsequences.length > 0 || filters.selectedImpacts.length > 0) {
            filteredVariant.molecularConsequences = variant.molecularConsequences?.filter(mc => {
                const matchesConsequence = filters.selectedConsequences.length === 0 ||
                    filters.selectedConsequences.includes(mc.featureConsequence || '');

                const matchesImpact = filters.selectedImpacts.length === 0 ||
                    filters.selectedImpacts.includes(mc.impact || '');

                return matchesConsequence && matchesImpact;
            });
        }

        return filteredVariant;
    });
}