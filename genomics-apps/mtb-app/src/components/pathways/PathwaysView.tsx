"use client";

import { useEffect, useMemo, useState } from 'react';
import { Variant } from '@/types/variants';
import { getVariantGene, getVariantShortLabel } from '@/lib/variantDisplay';
import { CooccurrenceMatch, CooccurrenceProfile } from '@/services/cooccurrenceService';
import { getDrawnGenes, getLoadedPathways, loadPathways, Pathway, PathwayNode } from '@/services/pathwayService';
import CooccurrenceBox from '../CooccurrenceBox';
import TxImplicationCell from '../TxImplicationCell';
import PathwayDiagram, { CooccurrenceConnector, PATHWAY_COLORS } from './PathwayDiagram';

interface PathwaysViewProps {
    selectedVariants: Variant[];
    cooccurrenceMatches: CooccurrenceMatch[];
    cooccurrencesByVariant: Record<string, CooccurrenceProfile[]>;
    onShowInResults: (variantId: string) => void;
}

interface ConnectorGroup extends CooccurrenceConnector {
    genes: string[];
    profiles: CooccurrenceProfile[];
}

const SOURCE_NOTE = 'Redrawn from Sanchez-Vega et al., Cell 2018 (Figure 2)';

export default function PathwaysView({
    selectedVariants,
    cooccurrenceMatches,
    cooccurrencesByVariant,
    onShowInResults,
}: PathwaysViewProps) {
    const [pathways, setPathways] = useState<Pathway[] | null>(getLoadedPathways);
    const [loadError, setLoadError] = useState<string | null>(null);
    const [pathwayId, setPathwayId] = useState<string | null>(null);
    const [hiddenVariantIds, setHiddenVariantIds] = useState<Set<string>>(new Set());
    const [selectedNodeId, setSelectedNodeId] = useState<string | null>(null);
    const [focusedConnectorKey, setFocusedConnectorKey] = useState<string | null>(null);

    useEffect(() => {
        loadPathways()
            .then(setPathways)
            .catch((error) => setLoadError(error instanceof Error ? error.message : 'Unable to load pathway diagrams'));
    }, []);

    const selectedGenes = useMemo(
        () => [...new Set(selectedVariants.map(getVariantGene).filter((gene): gene is string => Boolean(gene)))],
        [selectedVariants]
    );

    const rankedPathways = useMemo(() => (pathways ?? [])
        .map((pathway) => ({ pathway, hits: getDrawnGenes(pathway).filter((gene) => selectedGenes.includes(gene)) }))
        .sort((left, right) => right.hits.length - left.hits.length), [pathways, selectedGenes]);

    // Open on the top-ranked pathway
    useEffect(() => {
        if (!pathwayId && rankedPathways.length > 0) {
            setPathwayId(rankedPathways[0].pathway.id);
        }
    }, [pathwayId, rankedPathways]);

    const pathway = rankedPathways.find((entry) => entry.pathway.id === pathwayId)?.pathway ?? null;
    const diagramVariants = selectedVariants.filter((variant) => !variant.id || !hiddenVariantIds.has(variant.id));

    const variantsByNode = useMemo(() => {
        const byNode: Record<string, Variant[]> = {};
        pathway?.nodes.forEach((node) => {
            const variants = diagramVariants.filter((variant) => {
                const gene = getVariantGene(variant);
                return Boolean(gene && node.genes?.includes(gene));
            });
            if (variants.length > 0) {
                byNode[node.id] = variants;
            }
        });
        return byNode;
    }, [pathway, diagramVariants]);

    // One connector per pair of diagram boxes linked by matching co-occurrence profiles
    const connectorGroups = useMemo(() => {
        if (!pathway) {
            return [];
        }

        const nodeForVariant = (variantId: string) =>
            pathway.nodes.find((node) => (variantsByNode[node.id] ?? []).some((variant) => variant.id === variantId));
        const groups = new Map<string, ConnectorGroup>();

        cooccurrenceMatches.forEach((match) => {
            // Every component needs a selected, visible variant on this diagram
            const componentNodes = match.componentVariantIds.map((variantIds) =>
                variantIds.map(nodeForVariant).filter((node): node is PathwayNode => Boolean(node)));
            if (componentNodes.some((componentNodeList) => componentNodeList.length === 0)) {
                return;
            }
            const nodes = [...new Set(componentNodes.flat())];
            for (let i = 0; i < nodes.length; i += 1) {
                for (let j = i + 1; j < nodes.length; j += 1) {
                    const [first, second] = [nodes[i], nodes[j]].sort((a, b) => a.id.localeCompare(b.id));
                    const key = `${first.id}|${second.id}`;
                    const group = groups.get(key) ?? {
                        key,
                        fromNodeId: first.id,
                        toNodeId: second.id,
                        profileCount: 0,
                        genes: [first.label, second.label],
                        profiles: [],
                    };
                    if (!group.profiles.some((profile) => profile.id === match.profile.id)) {
                        group.profiles.push(match.profile);
                        group.profileCount = group.profiles.length;
                    }
                    groups.set(key, group);
                }
            }
        });

        return Array.from(groups.values());
    }, [pathway, variantsByNode, cooccurrenceMatches]);

    const focusedConnector = connectorGroups.find((group) => group.key === focusedConnectorKey) ?? null;
    const focusedGenes = useMemo(() => {
        if (!focusedConnector || !pathway) {
            return null;
        }
        const nodes = pathway.nodes.filter((node) => node.id === focusedConnector.fromNodeId || node.id === focusedConnector.toNodeId);
        return new Set(nodes.flatMap((node) => node.genes ?? []));
    }, [focusedConnector, pathway]);

    const selectPathway = (id: string) => {
        setPathwayId(id);
        setSelectedNodeId(null);
        setFocusedConnectorKey(null);
    };

    const toggleHidden = (variantId: string) => {
        setHiddenVariantIds((current) => {
            const next = new Set(current);
            if (next.has(variantId)) {
                next.delete(variantId);
            } else {
                next.add(variantId);
            }
            return next;
        });
    };

    const selectNode = (nodeId: string) => {
        setSelectedNodeId((current) => (current === nodeId ? null : nodeId));
        setFocusedConnectorKey(null);
    };

    const selectConnector = (key: string) => {
        setFocusedConnectorKey((current) => (current === key ? null : key));
        setSelectedNodeId(null);
    };

    const roleOf = (gene: string) => pathway?.genes.find((entry) => entry.symbol === gene)?.role ?? null;
    const notOnDiagram = pathway
        ? [...new Set(selectedVariants.map((variant) => getVariantGene(variant) ?? variant.range))]
            .filter((label) => !getDrawnGenes(pathway).includes(label))
        : [];
    const selectedNode = pathway?.nodes.find((node) => node.id === selectedNodeId) ?? null;
    const hits = rankedPathways.find((entry) => entry.pathway.id === pathwayId)?.hits ?? [];

    const panelClassName = 'min-w-0 rounded-xl border border-slate-300 bg-white p-4';
    const panelTitleClassName = 'mb-3 text-xs font-bold uppercase tracking-wide text-gray-500';

    if (loadError) {
        return <div className="rounded-xl border border-red-200 bg-red-50 p-4 text-sm text-red-700">{loadError}</div>;
    }

    if (!pathways || !pathway) {
        return <div className="p-4 text-sm text-gray-500">Loading pathway diagrams...</div>;
    }

    const renderSidePanel = () => {
        if (selectedNode) {
            const nodeVariants = variantsByNode[selectedNode.id] ?? [];
            const nodeGenes = selectedNode.genes ?? [];
            const roleLabels = [...new Set(nodeGenes.map(roleOf).filter(Boolean))]
                .map((role) => (role === 'OG' ? 'Oncogene' : 'Tumor suppressor gene'));
            const profiles = [...new Map(nodeVariants
                .flatMap((variant) => (variant.id ? cooccurrencesByVariant[variant.id] ?? [] : []))
                .map((profile) => [profile.id, profile])).values()];

            return (
                <div className={panelClassName}>
                    <button type="button" onClick={() => setSelectedNodeId(null)} className="text-sm font-semibold text-blue-600 hover:underline">
                        ← Back
                    </button>
                    <h3 className="mt-2 text-xl font-bold text-gray-900">{selectedNode.label}</h3>
                    <p className="mt-1 text-xs text-gray-500">
                        {nodeGenes.length > 1 ? `Genes: ${nodeGenes.join(', ')} · ` : ''}
                        {roleLabels.join(', ') || 'No role in Table S3'}
                    </p>
                    <div className="mt-4 border-t border-gray-200 pt-3">
                        <h4 className={panelTitleClassName}>Patient variants here</h4>
                        {nodeVariants.length === 0 && <p className="text-sm text-gray-500">No selected variants in this gene.</p>}
                        <div className="space-y-3">
                            {nodeVariants.map((variant) => (
                                <div key={variant.id ?? variant.variant} className="rounded-lg bg-gray-50 p-3">
                                    <div className="text-sm font-semibold text-gray-900">
                                        {getVariantGene(variant) ?? variant.range} {getVariantShortLabel(variant)}
                                    </div>
                                    <div className="mt-0.5 break-all font-mono text-xs text-gray-600">{variant.variant}</div>
                                    <div className="-mx-3">
                                        <TxImplicationCell implications={variant.txImplications} />
                                    </div>
                                    {variant.id && (
                                        <button
                                            type="button"
                                            onClick={() => onShowInResults(variant.id as string)}
                                            className="text-sm font-semibold text-blue-600 hover:underline"
                                        >
                                            Show in results table
                                        </button>
                                    )}
                                </div>
                            ))}
                        </div>
                    </div>
                    {profiles.length > 0 && (
                        <div className="mt-4 border-t border-gray-200 pt-3">
                            <CooccurrenceBox profiles={profiles} defaultOpen />
                        </div>
                    )}
                </div>
            );
        }

        return (
            <div className={panelClassName}>
                {connectorGroups.length > 0 && (
                    <div className="mb-5">
                        <h4 className={panelTitleClassName}>Potentially relevant co-occurring variants on this diagram</h4>
                        <div className="space-y-2">
                            {connectorGroups.map((group) => {
                                const isFocused = group.key === focusedConnectorKey;
                                return (
                                    <div
                                        key={group.key}
                                        className={`rounded-lg border p-3 ${isFocused ? 'border-violet-400 bg-violet-50' : 'border-gray-200 hover:border-violet-300'}`}
                                    >
                                        <button
                                            type="button"
                                            aria-pressed={isFocused}
                                            onClick={() => selectConnector(group.key)}
                                            className="w-full text-left"
                                        >
                                            <span className="text-xs font-bold uppercase tracking-wide text-violet-700">{group.genes.join(' + ')}</span>
                                            <span className="mt-1 block text-sm text-gray-800">{group.profiles.map((profile) => profile.name).join('; ')}</span>
                                        </button>
                                        {isFocused && (
                                            <div className="mt-3">
                                                <CooccurrenceBox profiles={group.profiles} defaultOpen />
                                            </div>
                                        )}
                                    </div>
                                );
                            })}
                        </div>
                    </div>
                )}
                <h4 className={panelTitleClassName}>About this diagram</h4>
                <p className="text-sm text-gray-600">
                    One of the 10 TCGA oncogenic signaling pathways, redrawn from Figure 2 of Sanchez-Vega et al., <i>Cell</i> 2018.
                    Gene roles are from the paper&apos;s Table S3. Click a gene to see this patient&apos;s variants in it.
                </p>
            </div>
        );
    };

    return (
        <div className="grid grid-cols-[260px_minmax(0,1fr)_340px] items-start gap-4">
            <aside className="flex min-w-0 flex-col gap-4">
                <div className={panelClassName}>
                    <h4 className={panelTitleClassName}>Pathways ranked by your genes</h4>
                    <div className="flex flex-col gap-2">
                        {rankedPathways.map(({ pathway: entry, hits: entryHits }) => {
                            const isActive = entry.id === pathwayId;
                            return (
                                <button
                                    key={entry.id}
                                    type="button"
                                    aria-pressed={isActive}
                                    onClick={() => selectPathway(entry.id)}
                                    className={`rounded-lg border px-3 py-2 text-left ${isActive ? 'border-blue-500 bg-blue-50 shadow-[inset_3px_0_0_0_rgb(37_99_235)]' : 'border-gray-200 hover:border-blue-400'}`}
                                >
                                    <span className="block text-sm font-semibold text-gray-900">{entry.name}</span>
                                    <span className="mt-1 block h-1.5 overflow-hidden rounded bg-gray-100">
                                        <span
                                            className="block h-full bg-blue-600"
                                            style={{ width: `${selectedGenes.length ? (entryHits.length / selectedGenes.length) * 100 : 0}%` }}
                                        />
                                    </span>
                                    <span className="mt-1 flex justify-between gap-2 text-xs text-gray-500">
                                        <span className="tabular-nums">{entryHits.length} of {selectedGenes.length} genes</span>
                                        <span className="truncate">{entryHits.join(', ') || '—'}</span>
                                    </span>
                                </button>
                            );
                        })}
                    </div>
                </div>
                <div className={panelClassName}>
                    <h4 className={panelTitleClassName}>Selected variants</h4>
                    <div className="flex flex-col gap-1">
                        {selectedVariants.map((variant) => {
                            const isHidden = Boolean(variant.id && hiddenVariantIds.has(variant.id));
                            const gene = getVariantGene(variant);
                            const role = gene ? roleOf(gene) : null;
                            return (
                                <label key={variant.id ?? variant.variant} className="flex cursor-pointer items-center gap-2 rounded px-1 py-1 hover:bg-gray-50">
                                    <input
                                        type="checkbox"
                                        checked={!isHidden}
                                        onChange={() => variant.id && toggleHidden(variant.id)}
                                        className="h-4 w-4 accent-blue-600"
                                    />
                                    <span className={`min-w-0 flex-1 text-sm ${isHidden ? 'text-gray-400 line-through' : 'text-gray-800'}`}>
                                        <b className="font-semibold">{gene ?? variant.range}</b> {getVariantShortLabel(variant)}
                                    </span>
                                    {role && (
                                        <span
                                            className="h-2.5 w-2.5 flex-shrink-0 rounded-full"
                                            style={{ background: role === 'OG' ? PATHWAY_COLORS.oncogeneStroke : PATHWAY_COLORS.tsgStroke }}
                                            title={role === 'OG' ? 'Oncogene' : 'Tumor suppressor gene'}
                                        />
                                    )}
                                </label>
                            );
                        })}
                    </div>
                    <p className="mt-2 text-xs text-gray-500">Untick a variant to hide it from the diagram without changing the selection.</p>
                </div>
                <div className={panelClassName}>
                    <h4 className={panelTitleClassName}>Not on this diagram</h4>
                    {notOnDiagram.length > 0 ? (
                        <div className="flex flex-wrap gap-1.5">
                            {notOnDiagram.map((label) => (
                                <span key={label} className="rounded-md border border-gray-200 bg-gray-50 px-2 py-0.5 text-xs font-semibold text-gray-700">{label}</span>
                            ))}
                        </div>
                    ) : (
                        <p className="text-sm text-gray-500">All selected genes appear on this diagram.</p>
                    )}
                </div>
            </aside>

            <section className="min-w-0 overflow-hidden rounded-xl border border-slate-300 bg-white">
                <div className="border-b border-gray-200 px-4 py-3">
                    <h3 className="text-lg font-bold text-gray-900">{pathway.name} pathway</h3>
                    <p className="text-sm text-gray-500">
                        {hits.length} of {selectedGenes.length} selected genes on this diagram · {SOURCE_NOTE}
                    </p>
                </div>
                <div className="flex flex-wrap items-center gap-x-5 gap-y-1.5 border-b border-gray-200 bg-gray-50 px-4 py-2 text-xs text-gray-600">
                    <LegendSwatch fill={PATHWAY_COLORS.oncogeneFill} stroke={PATHWAY_COLORS.oncogeneStroke} label="Oncogene" />
                    <LegendSwatch fill={PATHWAY_COLORS.tsgFill} stroke={PATHWAY_COLORS.tsgStroke} label="Tumor suppressor gene" />
                    <LegendSwatch fill="#ffffff" stroke={PATHWAY_COLORS.variantStroke} strokeWidth={3} label="Patient's selected variant" />
                    <span className="inline-flex items-center gap-1.5">
                        <svg width="28" height="12" aria-hidden="true">
                            <path d="M2,10 Q14,-2 26,10" fill="none" stroke={PATHWAY_COLORS.cooccurrence} strokeWidth={2.2} strokeDasharray="5 4" />
                        </svg>
                        Co-occurring variants
                    </span>
                    <span>→ activation &nbsp; ⊣ inhibition &nbsp; — part of complex</span>
                </div>
                <div className="overflow-x-auto p-2">
                    <PathwayDiagram
                        pathway={pathway}
                        variantsByNode={variantsByNode}
                        connectors={connectorGroups}
                        selectedNodeId={selectedNodeId}
                        focusedConnectorKey={focusedConnectorKey}
                        focusedGenes={focusedGenes}
                        onSelectNode={selectNode}
                        onSelectConnector={selectConnector}
                    />
                </div>
            </section>

            <aside className="min-w-0">{renderSidePanel()}</aside>
        </div>
    );
}

function LegendSwatch({ fill, stroke, label, strokeWidth = 2 }: { fill: string; stroke: string; label: string; strokeWidth?: number }) {
    return (
        <span className="inline-flex items-center gap-1.5">
            <span className="inline-block h-3.5 w-6 rounded-sm" style={{ background: fill, border: `${strokeWidth}px solid ${stroke}` }} />
            {label}
        </span>
    );
}
