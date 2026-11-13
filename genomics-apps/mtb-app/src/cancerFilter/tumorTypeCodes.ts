import Papa from "papaparse";

// The cancer-type codes that make an implication count for a tumor type, from MTB_TumorTypeImplicationCodes.csv.
// Rows in the "Any" category (Cancer, Solid Tumor) count for every tumor type.
export const PAN_CANCER_CATEGORY = "Any";

const DISEASE_ONTOLOGY = "https://disease-ontology.org";

type ImplicationCodeRow = {
    geneListCategory: string;
    system: string;
    code: string;
    display: string;
};

type CodingLike = { system?: string; code?: string };

// Tumor type (gene-list category) -> the system|code keys allowed for it
export type TumorTypeCodeIndex = Map<string, Set<string>>;

/**
 * A comparable key for a coding. Disease Ontology numbers are compared without leading zeros, since the genomics
 * server sends e.g. "DOID:50861" where CIViC and the reports have "DOID:0050861".
 */
export function codingKey(coding: CodingLike) {
    const code = (coding.code ?? "").trim();
    const normalized = coding.system === DISEASE_ONTOLOGY ? code.replace(/^(DOID:)0+(?=\d)/i, "$1").toUpperCase() : code;
    return `${coding.system ?? ""}|${normalized}`;
}

let codesPromise: Promise<TumorTypeCodeIndex> | null = null;

export function loadTumorTypeImplicationCodes() {
    if (!codesPromise) {
        codesPromise = new Promise<TumorTypeCodeIndex>((resolve, reject) => {
            Papa.parse("/data/MTB_TumorTypeImplicationCodes.csv", {
                header: true,
                download: true,
                skipEmptyLines: true,
                complete: (result) => {
                    const index: TumorTypeCodeIndex = new Map();
                    (result.data as ImplicationCodeRow[])
                        .filter((row) => row.geneListCategory && row.system && row.code)
                        .forEach((row) => {
                            const category = row.geneListCategory.trim();
                            const keys = index.get(category) ?? new Set<string>();
                            keys.add(codingKey({ system: row.system.trim(), code: row.code }));
                            index.set(category, keys);
                        });
                    resolve(index);
                },
                error: (error) => {
                    codesPromise = null;
                    reject(error);
                },
            });
        });
    }
    return codesPromise;
}

/**
 * True when any of the implication's cancer-type codings is allowed for the tumor type, or for every tumor type.
 * An implication with no codings never matches.
 */
export function matchesTumorType(codings: CodingLike[] | undefined, tumorType: string, index: TumorTypeCodeIndex) {
    if (!tumorType || !codings?.length) {
        return false;
    }
    const allowed = [index.get(tumorType), index.get(PAN_CANCER_CATEGORY)];
    return codings.some((coding) => allowed.some((keys) => keys?.has(codingKey(coding))));
}
