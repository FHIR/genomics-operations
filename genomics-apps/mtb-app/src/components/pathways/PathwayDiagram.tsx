import { KeyboardEvent } from 'react';
import { Variant } from '@/types/variants';
import { getVariantShortLabel } from '@/lib/variantDisplay';
import { Pathway, PathwayGroup, PathwayNode } from '@/services/pathwayService';

export const PATHWAY_COLORS = {
    oncogeneFill: '#f6d2e2',
    oncogeneStroke: '#d6337a',
    tsgFill: '#d1def5',
    tsgStroke: '#2d6bd0',
    neutralFill: '#ffffff',
    neutralStroke: '#9aa6b6',
    variantStroke: '#111827',
    edge: '#6b7787',
    group: '#97a3b4',
    cooccurrence: '#7c3aed',
    focus: '#2563eb',
};

export interface CooccurrenceConnector {
    key: string;
    fromNodeId: string;
    toNodeId: string;
    profileCount: number;
}

interface PathwayDiagramProps {
    pathway: Pathway;
    variantsByNode: Record<string, Variant[]>;
    connectors: CooccurrenceConnector[];
    selectedNodeId: string | null;
    focusedConnectorKey: string | null;
    focusedGenes: Set<string> | null;
    onSelectNode: (nodeId: string) => void;
    onSelectConnector: (key: string) => void;
}

type Box = PathwayNode | PathwayGroup;

const isGroup = (box: Box): box is PathwayGroup => 'members' in box;

const rectOf = (box: Box) => isGroup(box)
    ? [box.x, box.y, box.x + box.width, box.y + box.height]
    : [box.x - box.width / 2, box.y - box.height / 2, box.x + box.width / 2, box.y + box.height / 2];

const centerOf = (box: Box) => {
    const [left, top, right, bottom] = rectOf(box);
    return [(left + right) / 2, (top + bottom) / 2];
};

// Point where the line from the box center toward `target` leaves the box, plus a small gap
const clipTo = (box: Box, target: number[], gap: number) => {
    const [left, top, right, bottom] = rectOf(box);
    const [cx, cy] = centerOf(box);
    const dx = target[0] - cx;
    const dy = target[1] - cy;
    const t = Math.min(((right - left) / 2) / Math.abs(dx || 1e-9), ((bottom - top) / 2) / Math.abs(dy || 1e-9));
    const length = Math.hypot(dx, dy) || 1;
    return [cx + dx * t + (dx / length) * gap, cy + dy * t + (dy / length) * gap];
};

const activateOnKey = (event: KeyboardEvent, action: () => void) => {
    if (event.key === 'Enter' || event.key === ' ') {
        event.preventDefault();
        action();
    }
};

export default function PathwayDiagram({
    pathway,
    variantsByNode,
    connectors,
    selectedNodeId,
    focusedConnectorKey,
    focusedGenes,
    onSelectNode,
    onSelectConnector,
}: PathwayDiagramProps) {
    const boxes: Record<string, Box> = Object.fromEntries(
        [...pathway.nodes, ...pathway.groups].map((box) => [box.id, box])
    );
    const roles = Object.fromEntries(pathway.genes.map((gene) => [gene.symbol, gene.role]));
    const nodesById = Object.fromEntries(pathway.nodes.map((node) => [node.id, node]));
    const markerPrefix = `pathway-${pathway.id}`;

    return (
        <svg
            viewBox={`0 0 ${pathway.size.width} ${pathway.size.height}`}
            className="block h-auto w-full"
            style={{ minWidth: Math.round(pathway.size.width * 0.75) }}
            role="img"
            aria-label={`${pathway.name} pathway diagram`}
        >
            <defs>
                <marker id={`${markerPrefix}-activates`} viewBox="0 0 10 10" refX="9" refY="5" markerWidth="9" markerHeight="9" markerUnits="userSpaceOnUse" orient="auto">
                    <path d="M0,0 L10,5 L0,10 z" fill={PATHWAY_COLORS.edge} />
                </marker>
                <marker id={`${markerPrefix}-inhibits`} viewBox="0 0 4 14" refX="3" refY="7" markerWidth="4" markerHeight="14" markerUnits="userSpaceOnUse" orient="auto">
                    <rect x="0" y="0" width="4" height="14" fill={PATHWAY_COLORS.edge} />
                </marker>
            </defs>

            {pathway.groups.map((group) => {
                const position = group.labelPosition ?? 'top';
                const labelX = position === 'bottom-right' ? group.x + group.width : group.x + (position === 'top' ? 10 : 0);
                const labelY = position === 'top' ? group.y + 15 : group.y + group.height + 15;

                return (
                    <g key={group.id}>
                        <rect x={group.x} y={group.y} width={group.width} height={group.height} rx={6} fill="none" stroke={PATHWAY_COLORS.group} strokeWidth={1.2} />
                        {group.label && (
                            <text x={labelX} y={labelY} textAnchor={position === 'bottom-right' ? 'end' : 'start'} fontSize={12} fill="#4b5563">
                                {group.label}
                            </text>
                        )}
                    </g>
                );
            })}

            {pathway.edges.map((edge, index) => {
                const from = boxes[edge.from];
                const to = boxes[edge.to];
                if (!from || !to) {
                    return null;
                }

                const [x1, y1] = clipTo(from, centerOf(to), 3);
                const [x2, y2] = clipTo(to, centerOf(from), edge.type === 'binds' ? 3 : 4);
                const marker = edge.type === 'binds' ? undefined : `url(#${markerPrefix}-${edge.type})`;

                return <line key={`${edge.from}-${edge.to}-${index}`} x1={x1} y1={y1} x2={x2} y2={y2} stroke={PATHWAY_COLORS.edge} strokeWidth={1.6} markerEnd={marker} />;
            })}

            {pathway.nodes.map((node) => {
                const [left, top] = rectOf(node);

                if (node.type !== 'gene') {
                    return (
                        <text key={node.id} x={node.x} y={node.y + 4.5} textAnchor="middle" fontSize={13} fontWeight={500} fill="#4b5563">
                            {node.label}
                        </text>
                    );
                }

                const role = (node.genes ?? []).map((gene) => roles[gene]).find(Boolean);
                const fill = role === 'OG' ? PATHWAY_COLORS.oncogeneFill : role === 'TSG' ? PATHWAY_COLORS.tsgFill : PATHWAY_COLORS.neutralFill;
                const roleStroke = role === 'OG' ? PATHWAY_COLORS.oncogeneStroke : role === 'TSG' ? PATHWAY_COLORS.tsgStroke : PATHWAY_COLORS.neutralStroke;
                const variants = variantsByNode[node.id] ?? [];
                const hasVariant = variants.length > 0;
                const isSelected = selectedNodeId === node.id;
                const isFocused = Boolean(focusedGenes && (node.genes ?? []).some((gene) => focusedGenes.has(gene)));
                const isDimmed = Boolean(focusedGenes) && !isFocused;
                const variantLabel = variants.map(getVariantShortLabel).join(', ');

                return (
                    <g
                        key={node.id}
                        role="button"
                        tabIndex={0}
                        aria-label={`${node.label}${hasVariant ? `, ${variantLabel}` : ''}`}
                        onClick={() => onSelectNode(node.id)}
                        onKeyDown={(event) => activateOnKey(event, () => onSelectNode(node.id))}
                        className="cursor-pointer focus:outline-none"
                        opacity={isDimmed ? 0.35 : 1}
                    >
                        {(isSelected || isFocused) && (
                            <rect x={left - 4} y={top - 4} width={node.width + 8} height={node.height + 8} rx={9} fill="none" stroke={PATHWAY_COLORS.focus} strokeWidth={2.5} strokeDasharray="6 4" />
                        )}
                        <rect
                            x={left}
                            y={top}
                            width={node.width}
                            height={node.height}
                            rx={6}
                            fill={fill}
                            stroke={hasVariant ? PATHWAY_COLORS.variantStroke : roleStroke}
                            strokeWidth={hasVariant ? 3 : 1.4}
                        />
                        {hasVariant ? (
                            <>
                                <text x={node.x} y={node.y - 3} textAnchor="middle" fontSize={13} fontWeight={700} fill="#111827">{node.label}</text>
                                <text x={node.x} y={node.y + 12} textAnchor="middle" fontSize={10} fontFamily="ui-monospace, SFMono-Regular, Menlo, monospace" fill="#111827">
                                    {variantLabel.length > 15 ? `${variantLabel.slice(0, 14)}…` : variantLabel}
                                </text>
                            </>
                        ) : (
                            <text x={node.x} y={node.y + 4.5} textAnchor="middle" fontSize={13} fontWeight={700} fill="#111827">{node.label}</text>
                        )}
                        {variants.length > 1 && (
                            <>
                                <circle cx={left + node.width} cy={top} r={9} fill="#111827" />
                                <text x={left + node.width} y={top + 4} textAnchor="middle" fontSize={11} fontWeight={700} fill="#ffffff">{variants.length}</text>
                            </>
                        )}
                    </g>
                );
            })}

            {connectors.map((connector) => {
                const a = nodesById[connector.fromNodeId];
                const b = nodesById[connector.toNodeId];
                if (!a || !b) {
                    return null;
                }

                // Arc above the two genes, kept inside the drawing area
                const dx = b.x - a.x;
                const dy = b.y - a.y;
                const length = Math.hypot(dx, dy) || 1;
                let px = dy / length;
                let py = -dx / length;
                if (py > 0) {
                    px = -px;
                    py = -py;
                }
                const top = (a.y - a.height / 2 + b.y - b.height / 2) / 2;
                let offset = Math.max(40, 0.35 * length);
                if (py < 0) {
                    offset = Math.min(offset, (2 * (top - 14)) / -py);
                }
                const cx = (a.x + b.x) / 2 + px * offset;
                const cy = top + py * offset;
                const midX = 0.25 * a.x + 0.5 * cx + 0.25 * b.x;
                const midY = 0.5 * top + 0.5 * cy;
                const path = `M${a.x},${a.y - a.height / 2} Q${cx},${cy} ${b.x},${b.y - b.height / 2}`;
                const isFocusedConnector = focusedConnectorKey === connector.key;

                return (
                    <g
                        key={connector.key}
                        role="button"
                        tabIndex={0}
                        aria-label={`Co-occurring variants: ${a.label} and ${b.label}`}
                        onClick={() => onSelectConnector(connector.key)}
                        onKeyDown={(event) => activateOnKey(event, () => onSelectConnector(connector.key))}
                        className="cursor-pointer focus:outline-none"
                    >
                        <path d={path} fill="none" stroke="transparent" strokeWidth={16} />
                        <path d={path} fill="none" stroke={PATHWAY_COLORS.cooccurrence} strokeWidth={isFocusedConnector ? 4 : 2.4} strokeDasharray="7 5" pointerEvents="none" />
                        <circle cx={midX} cy={midY} r={10} fill={PATHWAY_COLORS.cooccurrence} />
                        <text x={midX} y={midY + 4} textAnchor="middle" fontSize={11} fontWeight={700} fill="#ffffff" pointerEvents="none">{connector.profileCount}</text>
                    </g>
                );
            })}
        </svg>
    );
}
