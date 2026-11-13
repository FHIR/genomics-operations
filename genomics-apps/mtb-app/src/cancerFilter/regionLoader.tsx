import { useEffect, useState } from "react";
import { KBRow, PresetSelection, loadGeneListKbRows } from "./geneListKb";

interface RegionLoaderProps {
    presetSelection: PresetSelection | null;
    requestId: number;
    onRegionsLoaded: (regions: string[]) => void;
}

export default function RegionLoader({ presetSelection, requestId, onRegionsLoaded }: RegionLoaderProps) {
    const [kbRows, setKbRows] = useState<KBRow[]>([]);

    useEffect(() => {
        loadGeneListKbRows().then(setKbRows);
    }, []);

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
