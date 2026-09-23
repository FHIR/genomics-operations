import Papa from "papaparse";

export const GLOBAL_GENE_LIST_CATEGORY = "Any";

export type KBRow = {
    "Cancer category": string;
    "Phenotype (for filtering implications)": string;
    "Label": string;
    "Region": string;
};

export type PresetSelection = {
    category: string;
    label: string;
};

function normalizeRow(row: KBRow): KBRow {
    return {
        "Cancer category": row["Cancer category"]?.trim() ?? "",
        "Phenotype (for filtering implications)": row["Phenotype (for filtering implications)"]?.trim() ?? "",
        "Label": row["Label"]?.trim() ?? "",
        "Region": row["Region"]?.trim() ?? "",
    };
}

export function loadGeneListKbRows() {
    return new Promise<KBRow[]>((resolve, reject) => {
        Papa.parse("/data/MTB_KB_GeneLists.csv", {
            header: true,
            download: true,
            complete: (result) => {
                const rows = (result.data as KBRow[])
                    .map(normalizeRow)
                    .filter((row) => row["Cancer category"] || row["Label"] || row["Region"]);

                resolve(rows);
            },
            error: (error) => {
                reject(error);
            },
        });
    });
}

export function getAvailableCancerTypes(rows: KBRow[]) {
    return Array.from(
        new Set(
            rows
                .map((row) => row["Cancer category"])
                .filter((category): category is string => Boolean(category) && category !== GLOBAL_GENE_LIST_CATEGORY)
        )
    ).sort((left, right) => left.localeCompare(right));
}

export function getGlobalGeneListLabels(rows: KBRow[]) {
    return Array.from(
        new Set(
            rows
                .filter((row) => row["Cancer category"] === GLOBAL_GENE_LIST_CATEGORY)
                .map((row) => row["Label"])
                .filter((label): label is string => Boolean(label))
        )
    ).sort((left, right) => left.localeCompare(right));
}

export function getCancerSpecificPresetLabels(rows: KBRow[], cancerType: string) {
    if (!cancerType) {
        return [];
    }

    const globalLabels = new Set(getGlobalGeneListLabels(rows));

    return Array.from(
        new Set(
            rows
                .filter((row) => row["Cancer category"] === cancerType)
                .map((row) => row["Label"])
                .filter((label): label is string => Boolean(label) && !globalLabels.has(label))
        )
    ).sort((left, right) => left.localeCompare(right));
}
