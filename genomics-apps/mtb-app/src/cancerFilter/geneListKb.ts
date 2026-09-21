import Papa from "papaparse";

export type KBRow = {
    "Cancer category": string;
    "Phenotype (for filtering implications)": string;
    "Label": string;
    "Region": string;
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
