import { useEffect, useState } from "react";
import { KBRow, PresetSelection, loadGeneListKbRows } from "./geneListKb";

interface RegionLoaderProps {
    presetSelection: PresetSelection | null;
    requestId: number;
    onRegionsLoaded: (regions: string[]) => void;
    onPhenotypesLoaded?: (phenotypes: Set<string>) => void;
}



export default function RegionLoader({ presetSelection, requestId, onRegionsLoaded, onPhenotypesLoaded }: RegionLoaderProps) {
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
        if (!presetSelection || requestId === 0 || kbRows.length === 0) return;

        // Load predefined gene lists for the selected preset category and label.
        const matchedRegions = kbRows
            .filter((row) => row["Cancer category"] === presetSelection.category && row["Label"] === presetSelection.label)
            .map((row) => row["Region"]);

        onRegionsLoaded(matchedRegions);
    }, [presetSelection, requestId, kbRows, onRegionsLoaded]);

    return null; // it only loads data
}
