// Helper function to check phenotype match for individual implications
export function isPhenotypeMatchForImplication(
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

export function isOtherTumorType(
  phenotype: string | undefined,
  selectedCancerType: string,
  mtbKbPhenotypes: Set<string>
): boolean {
  if (!phenotype || !mtbKbPhenotypes.has(phenotype)) return false;
  
  // Check if it's NOT the selected cancer type
  const selectedCancerPhenotypes = new Set<string>();
  if (mtbKbPhenotypes.size > 0) {
    mtbKbPhenotypes.forEach(p => {
      if (selectedCancerType === "Breast carcinoma" && p.toLowerCase().includes("breast")) {
        selectedCancerPhenotypes.add(p);
      } else if (selectedCancerType === "NSCLC" && p.toLowerCase().includes("lung")) {
        selectedCancerPhenotypes.add(p);
      }
    });
  }
  
  return !selectedCancerPhenotypes.has(phenotype);
}

// Helper function to check if evidence level is A-level
export function isALevelEvidence(evidenceLevel: string | undefined): boolean {
  if (!evidenceLevel) return false;
  
  return evidenceLevel.startsWith("A ") || 
         evidenceLevel === "A" || 
         evidenceLevel.includes("Validated association");
}