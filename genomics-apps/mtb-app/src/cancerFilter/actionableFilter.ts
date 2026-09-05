import { Variant } from '@/types/variants';

// Helper function to check phenotype match for individual implications
function isPhenotypeMatchForImplication(
  phenotype: string | undefined,
  selectedCancerType: string
): boolean {
  if (!phenotype || !selectedCancerType) return false;

  // Split by semicolon and check each phenotype
  const phenotypes = phenotype.split(';').map(p => p.trim());
  return phenotypes.some(p => {
    const phenotypeLower = p.toLowerCase();
    const cancerTypeLower = selectedCancerType.toLowerCase();

    // For NSCLC, also check for "lung" and "non-small cell" patterns
    if (selectedCancerType === "NSCLC") {
      return phenotypeLower.includes("lung") && phenotypeLower.includes("non-small");
    }

    // For Breast carcinoma, check for "breast" patterns
    if (selectedCancerType === "Breast carcinoma") {
      return phenotypeLower.includes("breast");
    }

    // For other cancer types, use direct matching
    return phenotypeLower.includes(cancerTypeLower);
  });
}

export function filterActionableThisTumor(
  items: Variant[],
  selectedCancerType: string,
  mtbKbPhenotypes?: Set<string>
): Variant[] {
  console.log(`filterActionableThisTumor called with:`, {
    selectedCancerType,
    mtbKbPhenotypesSize: mtbKbPhenotypes?.size,
    itemsCount: items.length
  });

  return items.map(item => {
    console.log(`\n=== Filtering item: ${item.variant} ===`);

    // Filter TX implications to only include those with valid evidence and matching phenotype
    const filteredTxImplications = item.txImplications?.filter(tx => {
      // Check evidence level
      const hasValidEvidence = tx.evidenceLevel?.startsWith("A ") || tx.evidenceLevel === "A" || tx.evidenceLevel?.includes("Validated association");
      console.log(`TX Evidence check: ${tx.evidenceLevel} -> ${hasValidEvidence}`);

      if (!hasValidEvidence) return false;

      // Check phenotype match
      let phenotypeMatch: boolean;
      if (mtbKbPhenotypes && mtbKbPhenotypes.size > 0) {
        // Check if phenotype is in MTB KB AND matches the selected cancer type
        const phenotypeInKB = tx.phenotypicContext && mtbKbPhenotypes.has(tx.phenotypicContext);
        const phenotypeMatchesCancerType = isPhenotypeMatchForImplication(tx.phenotypicContext, selectedCancerType);
        phenotypeMatch = !!(phenotypeInKB && phenotypeMatchesCancerType);
        console.log(`TX phenotype in KB: ${tx.phenotypicContext} -> ${phenotypeInKB}`);
        console.log(`TX phenotype matches cancer type: ${tx.phenotypicContext} -> ${phenotypeMatchesCancerType}`);
        console.log(`TX final match: ${phenotypeMatch}`);
      } else {
        phenotypeMatch = isPhenotypeMatchForImplication(tx.phenotypicContext, selectedCancerType);
        console.log(`TX fallback phenotype match: ${tx.phenotypicContext} -> ${phenotypeMatch}`);
      }

      return phenotypeMatch;
    }) || [];

    // Filter DX implications to only include those with valid evidence and matching phenotype
    const filteredDxImplications = item.dxImplications?.filter(dx => {
      // Check evidence level
      const hasValidEvidence = dx.evidenceLevel?.startsWith("A ") || dx.evidenceLevel === "A" || dx.evidenceLevel?.includes("Validated association");
      console.log(`DX Evidence check: ${dx.evidenceLevel} -> ${hasValidEvidence}`);

      if (!hasValidEvidence) return false;

      // Check phenotype match
      let phenotypeMatch: boolean;
      if (mtbKbPhenotypes && mtbKbPhenotypes.size > 0) {
        // Check if phenotype is in MTB KB AND matches the selected cancer type
        const phenotypeInKB = dx.predictedPhenotype && mtbKbPhenotypes.has(dx.predictedPhenotype);
        const phenotypeMatchesCancerType = isPhenotypeMatchForImplication(dx.predictedPhenotype, selectedCancerType);
        phenotypeMatch = !!(phenotypeInKB && phenotypeMatchesCancerType);
        console.log(`DX phenotype in KB: ${dx.predictedPhenotype} -> ${phenotypeInKB}`);
        console.log(`DX phenotype matches cancer type: ${dx.predictedPhenotype} -> ${phenotypeMatchesCancerType}`);
        console.log(`DX final match: ${phenotypeMatch}`);
      } else {
        phenotypeMatch = isPhenotypeMatchForImplication(dx.predictedPhenotype, selectedCancerType);
        console.log(`DX fallback phenotype match: ${dx.predictedPhenotype} -> ${phenotypeMatch}`);
      }

      return phenotypeMatch;
    }) || [];

    // Only include the variant if it has at least one matching implication
    if (filteredTxImplications.length > 0 || filteredDxImplications.length > 0) {
      return {
        ...item,
        txImplications: filteredTxImplications,
        dxImplications: filteredDxImplications
      };
    }

    return null;
  }).filter((item): item is Variant => item !== null);
}

export function filterActionableOtherTumor(
  items: Variant[],
  selectedCancerType: string,
  mtbKbPhenotypes: Set<string>
): Variant[] {
  console.log(`filterActionableOtherTumor called with:`, {
    selectedCancerType,
    mtbKbPhenotypesSize: mtbKbPhenotypes.size,
    itemsCount: items.length
  });

  // Get phenotypes for the selected cancer type to exclude them
  const selectedCancerPhenotypes = new Set<string>();
  if (mtbKbPhenotypes.size > 0) {
    // We need to get the specific phenotypes for the selected cancer type
    // For now, we'll use a simple approach based on the cancer type
    mtbKbPhenotypes.forEach(phenotype => {
      if (selectedCancerType === "Breast carcinoma" &&
          (phenotype.toLowerCase().includes("breast"))) {
        selectedCancerPhenotypes.add(phenotype);
      } else if (selectedCancerType === "NSCLC" &&
                 (phenotype.toLowerCase().includes("lung"))) {
        selectedCancerPhenotypes.add(phenotype);
      }
    });
  }

  console.log('Selected cancer phenotypes to exclude:', Array.from(selectedCancerPhenotypes));

  return items.map(item => {
    console.log(`\n=== Filtering item: ${item.variant} ===`);

    // Filter TX implications to only include those with valid evidence and phenotypes from OTHER tumor types
    const filteredTxImplications = item.txImplications?.filter(tx => {
      // Check evidence level
      const hasValidEvidence = tx.evidenceLevel?.startsWith("A ") || tx.evidenceLevel === "A" || tx.evidenceLevel?.includes("Validated association");
      console.log(`TX Evidence check: ${tx.evidenceLevel} -> ${hasValidEvidence}`);

      if (!hasValidEvidence) return false;

      // Check if phenotype is in MTB KB but NOT for the selected cancer type
      const phenotypeInKB = tx.phenotypicContext && mtbKbPhenotypes.has(tx.phenotypicContext);
      const phenotypeNotForSelectedCancer = tx.phenotypicContext && !selectedCancerPhenotypes.has(tx.phenotypicContext);
      const isOtherTumorType = phenotypeInKB && phenotypeNotForSelectedCancer;

      console.log(`TX phenotype in KB: ${tx.phenotypicContext} -> ${phenotypeInKB}`);
      console.log(`TX phenotype not for selected cancer: ${tx.phenotypicContext} -> ${phenotypeNotForSelectedCancer}`);
      console.log(`TX is other tumor type: ${isOtherTumorType}`);

      return isOtherTumorType;
    }) || [];

    // Filter DX implications to only include those with valid evidence and phenotypes from OTHER tumor types
    const filteredDxImplications = item.dxImplications?.filter(dx => {
      // Check evidence level
      const hasValidEvidence = dx.evidenceLevel?.startsWith("A ") || dx.evidenceLevel === "A" || dx.evidenceLevel?.includes("Validated association");
      console.log(`DX Evidence check: ${dx.evidenceLevel} -> ${hasValidEvidence}`);

      if (!hasValidEvidence) return false;

      // Check if phenotype is in MTB KB but NOT for the selected cancer type
      const phenotypeInKB = dx.predictedPhenotype && mtbKbPhenotypes.has(dx.predictedPhenotype);
      const phenotypeNotForSelectedCancer = dx.predictedPhenotype && !selectedCancerPhenotypes.has(dx.predictedPhenotype);
      const isOtherTumorType = phenotypeInKB && phenotypeNotForSelectedCancer;

      console.log(`DX phenotype in KB: ${dx.predictedPhenotype} -> ${phenotypeInKB}`);
      console.log(`DX phenotype not for selected cancer: ${dx.predictedPhenotype} -> ${phenotypeNotForSelectedCancer}`);
      console.log(`DX is other tumor type: ${isOtherTumorType}`);

      return isOtherTumorType;
    }) || [];

    // Only include the variant if it has at least one matching implication
    if (filteredTxImplications.length > 0 || filteredDxImplications.length > 0) {
      return {
        ...item,
        txImplications: filteredTxImplications,
        dxImplications: filteredDxImplications
      };
    }

    return null;
  }).filter((item): item is Variant => item !== null);
}

export function filterPossiblyActionable(
  items: Variant[]
): Variant[] {
  console.log(`filterPossiblyActionable called with:`, {
    itemsCount: items.length
  });

  return items.map(item => {
    console.log(`\n=== Filtering item: ${item.variant} ===`);

    // Filter TX implications to only include those WITHOUT evidence level A
    const filteredTxImplications = item.txImplications?.filter(tx => {
      // Check if evidence level is NOT A level
      const hasValidEvidence = tx.evidenceLevel?.startsWith("A ") || tx.evidenceLevel === "A" || tx.evidenceLevel?.includes("Validated association");
      const isPossiblyActionable = !hasValidEvidence;
      console.log(`TX Evidence check: ${tx.evidenceLevel} -> A level: ${hasValidEvidence}, possibly actionable: ${isPossiblyActionable}`);

      return isPossiblyActionable;
    }) || [];

    // Filter DX implications to only include those WITHOUT evidence level A
    const filteredDxImplications = item.dxImplications?.filter(dx => {
      // Check if evidence level is NOT A level
      const hasValidEvidence = dx.evidenceLevel?.startsWith("A ") || dx.evidenceLevel === "A" || dx.evidenceLevel?.includes("Validated association");
      const isPossiblyActionable = !hasValidEvidence;
      console.log(`DX Evidence check: ${dx.evidenceLevel} -> A level: ${hasValidEvidence}, possibly actionable: ${isPossiblyActionable}`);

      return isPossiblyActionable;
    }) || [];

    // Only include the variant if it has at least one matching implication
    if (filteredTxImplications.length > 0 || filteredDxImplications.length > 0) {
      return {
        ...item,
        txImplications: filteredTxImplications,
        dxImplications: filteredDxImplications
      };
    }

    return null;
  }).filter((item): item is Variant => item !== null);
}