import { useEffect, useState } from "react";
import { KBRow, loadGeneListKbRows } from "./geneListKb";

interface RegionLoaderProps {
    cancerType: string;
    label: string;
    requestId: number;
    onRegionsLoaded: (regions: string[]) => void;
    onPhenotypesLoaded?: (phenotypes: Set<string>) => void;
}



export default function RegionLoader({ cancerType, label, requestId, onRegionsLoaded, onPhenotypesLoaded }: RegionLoaderProps) {
    const [kbRows, setKbRows] = useState<KBRow[]>([]);

    useEffect(() => {
        loadGeneListKbRows().then((rows) => {
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
        });
    }, [onPhenotypesLoaded]);

    useEffect(() => {
        if (!cancerType || !label || requestId === 0 || kbRows.length === 0) return;

        // Load predefined gene lists for the selected cancer type and preset label.
        const matchedRegions = kbRows
            .filter((row) => row["Cancer category"] === cancerType && row["Label"] === label)
            .map((row) => row["Region"]);

        onRegionsLoaded(matchedRegions);
    }, [cancerType, label, requestId, kbRows, onRegionsLoaded]);

    return null; // it only loads data
}
