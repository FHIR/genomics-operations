// A-level evidence: the level text starts with "A", or is CIViC's "Validated association"
export function isALevelEvidence(evidenceLevel: string | undefined): boolean {
  if (!evidenceLevel) return false;

  return evidenceLevel.startsWith("A ") ||
    evidenceLevel === "A" ||
    evidenceLevel.includes("Validated association");
}
