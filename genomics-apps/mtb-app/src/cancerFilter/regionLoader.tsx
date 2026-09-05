import Papa from "papaparse";
import { useEffect, useState } from "react";

type KBRow = {
    "Cancer category": string;
    "Phenotype (for filtering implications)": string;
    "Label": string;
    "Region": string;
};

interface RegionLoaderProps {
    cancerType: string;
    label: string;
    onRegionsLoaded: (regions: string[]) => void;
    onPhenotypesLoaded?: (phenotypes: Set<string>) => void;
}



export default function RegionLoader({ cancerType, label, onRegionsLoaded, onPhenotypesLoaded }: RegionLoaderProps) {
    const [kbRows, setKbRows] = useState<KBRow[]>([]);

    useEffect(() => {
        Papa.parse("/data/MTB_KB.csv", {
            header: true,
            download: true,
            complete: (result) => {
                const rows = result.data as KBRow[];
                setKbRows(rows);

                // Extract phenotypes if callback is provided
                if (onPhenotypesLoaded) {
                    const phenotypes = new Set<string>();
                    rows.forEach(row => {
                        const phenotypeColumn = row["Phenotype (for filtering implications)"];
                        if (phenotypeColumn) {
                            const individualPhenotypes = phenotypeColumn.split(';').map(p => p.trim());
                            individualPhenotypes.forEach(phenotype => {
                                if (phenotype) {
                                    phenotypes.add(phenotype);
                                }
                            });
                        }
                    });
                    onPhenotypesLoaded(phenotypes);
                }
            },
        });
    }, [onPhenotypesLoaded]);

    useEffect(() => {
        if (!cancerType || !label || kbRows.length === 0) return;

        // Use cancer type directly from CSV data
        const matchedRegions = kbRows
            .filter((row) => row["Cancer category"] === cancerType && row["Label"] === label)
            .map((row) => row["Region"]);

        onRegionsLoaded(matchedRegions);
    }, [cancerType, label, kbRows, onRegionsLoaded]);

    return null; // it only loads data
}
