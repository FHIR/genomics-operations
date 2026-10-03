// PATHWAY DIAGRAMS (TCGA oncogenic signaling pathways stored in public/data/pathways)
const PATHWAY_DATA_PATH = '/data/pathways';

export interface PathwayNode {
    id: string;
    label: string;
    type: 'gene' | 'process';
    x: number; // center
    y: number;
    width: number;
    height: number;
    genes?: string[];
}

export interface PathwayGroup {
    id: string;
    label: string;
    x: number; // top-left corner
    y: number;
    width: number;
    height: number;
    members: string[];
    labelPosition?: 'top' | 'bottom-left' | 'bottom-right';
}

export interface PathwayEdge {
    from: string;
    to: string;
    type: 'activates' | 'inhibits' | 'binds';
}

export interface PathwayGene {
    symbol: string;
    role: 'OG' | 'TSG' | null;
    inDiagram: boolean;
}

export interface Pathway {
    id: string;
    name: string;
    size: { width: number; height: number };
    nodes: PathwayNode[];
    groups: PathwayGroup[];
    edges: PathwayEdge[];
    genes: PathwayGene[];
}

interface PathwayIndex {
    pathways: { id: string; name: string; file: string }[];
}

let pathwaysPromise: Promise<Pathway[]> | null = null;
let loadedPathways: Pathway[] | null = null;

// Pathways already loaded in this page view, so the Pathways view can render them immediately
export const getLoadedPathways = () => loadedPathways;

// Loads all pathway diagrams once per page view
export function loadPathways() {
    if (!pathwaysPromise) {
        pathwaysPromise = (async () => {
            const indexResponse = await fetch(`${PATHWAY_DATA_PATH}/index.json`);
            if (!indexResponse.ok) {
                throw new Error(`Unable to load pathway index (${indexResponse.status})`);
            }

            const index = await indexResponse.json() as PathwayIndex;
            return Promise.all(index.pathways.map(async (entry) => {
                const response = await fetch(`${PATHWAY_DATA_PATH}/${entry.file}`);
                if (!response.ok) {
                    throw new Error(`Unable to load pathway ${entry.id} (${response.status})`);
                }

                return response.json() as Promise<Pathway>;
            }));
        })().then((pathways) => {
            loadedPathways = pathways;
            return pathways;
        }).catch((error) => {
            pathwaysPromise = null;
            throw error;
        });
    }

    return pathwaysPromise;
}

// Genes drawn in the diagram; used for ranking pathways against the selected variants
export const getDrawnGenes = (pathway: Pathway) =>
    pathway.genes.filter((gene) => gene.inDiagram).map((gene) => gene.symbol);
