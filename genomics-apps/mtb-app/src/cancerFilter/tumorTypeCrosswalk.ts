import Papa from "papaparse";
import { Coding } from "@/types/fhir";

// Generated from the VSAC value sets by terminology/build_crosswalk.py
export type CrosswalkRow = {
    system: string;
    code: string;
    display: string;
    predictorTumorType: string;
    geneListCategory: string;
    valueSetOid: string;
};

export type TumorTypeMapping = {
    // Sent to the oncogenicity predictor; undefined when the predictor has no matching type
    predictorTumorType?: string;
    // Gene-list "Cancer category" for presets and implication filtering
    geneListCategory?: string;
};

// One entry of the banner's Dx dropdown
export type CaseDx = TumorTypeMapping & {
    key: string;
    label: string;
    source: "patient" | "generic";
    // The primary cancer Condition, for entries from the patient's record
    conditionId?: string;
};

let crosswalkPromise: Promise<CrosswalkRow[]> | null = null;

export function loadTumorTypeCrosswalk() {
    if (!crosswalkPromise) {
        crosswalkPromise = new Promise<CrosswalkRow[]>((resolve, reject) => {
            Papa.parse("/data/MTB_TumorTypeCrosswalk.csv", {
                header: true,
                download: true,
                skipEmptyLines: true,
                complete: (result) => resolve((result.data as CrosswalkRow[]).filter((row) => row.system && row.code)),
                error: (error) => {
                    crosswalkPromise = null;
                    reject(error);
                },
            });
        });
    }
    return crosswalkPromise;
}

/**
 * Exact code match of a condition's codings against the crosswalk. Codes under a listed code don't count.
 */
export function mapCodingsToTumorType(codings: Coding[] | undefined, rows: CrosswalkRow[]): TumorTypeMapping {
    for (const coding of codings ?? []) {
        const row = rows.find((candidate) => candidate.system === coding.system && candidate.code === coding.code);
        if (row) {
            return {
                predictorTumorType: row.predictorTumorType || undefined,
                geneListCategory: row.geneListCategory || undefined,
            };
        }
    }
    return {};
}

/**
 * The generic tumor types the crosswalk knows, labeled by gene-list category, alphabetical.
 */
export function getGenericDxOptions(rows: CrosswalkRow[]): CaseDx[] {
    const byCategory = new Map<string, CaseDx>();
    for (const row of rows) {
        if (row.geneListCategory && !byCategory.has(row.geneListCategory)) {
            byCategory.set(row.geneListCategory, {
                key: `generic:${row.geneListCategory}`,
                label: row.geneListCategory,
                source: "generic",
                predictorTumorType: row.predictorTumorType || undefined,
                geneListCategory: row.geneListCategory,
            });
        }
    }
    return [...byCategory.values()].sort((left, right) => left.label.localeCompare(right.label));
}
