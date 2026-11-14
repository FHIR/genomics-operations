import { StructuralDetails, Variant } from '@/types/variants';

const GENE_SYMBOL_PATTERN = /^[A-Z0-9-]+$/;

const AMINO_ACIDS: Record<string, string> = {
    Ala: 'A', Arg: 'R', Asn: 'N', Asp: 'D', Cys: 'C', Gln: 'Q', Glu: 'E', Gly: 'G', His: 'H', Ile: 'I',
    Leu: 'L', Lys: 'K', Met: 'M', Phe: 'F', Pro: 'P', Ser: 'S', Thr: 'T', Trp: 'W', Tyr: 'Y', Val: 'V',
    Ter: '*', Sec: 'U', Pyl: 'O',
};

export const getSimpleVariantLabel = (variantString: string) => {
    const parts = variantString.split(':');

    if (parts.length !== 4) {
        return 'Simple';
    }

    const [, , deleted = '', inserted = ''] = parts;

    if (deleted.length === inserted.length) {
        if (deleted.length === 1) {
            return 'SNV';
        }

        if (deleted.length > 1) {
            return 'MNV';
        }
    }

    if (deleted.length !== inserted.length) {
        return 'InDel';
    }

    return 'Simple';
};

// The gene a variant belongs to, when the search term was a gene symbol
export const getVariantGene = (variant: Variant) =>
    GENE_SYMBOL_PATTERN.test(variant.range) ? variant.range : undefined;

// Structural variants and fusions, described from their recorded fields (the report card and the Explorer share these)

// "copy_number_variation" becomes "Copy number variation"
export const getChangeTypeLabel = (details: StructuralDetails) => {
    const text = details.changeType?.replace(/_/g, ' ').trim();
    return text ? text.charAt(0).toUpperCase() + text.slice(1) : 'Structural variant';
};

// "1 copy", "4 copies"; absent when the copy number isn't recorded
export const getCopiesText = (details: StructuralDetails) =>
    details.copyNumber === undefined ? undefined : `${details.copyNumber} ${details.copyNumber === 1 ? 'copy' : 'copies'}`;

// "NC_000006.11:117609462-117747018", or just the reference sequence when there's no range
export const getStructuralLocation = (details: StructuralDetails) => {
    if (!details.referenceSequence) {
        return undefined;
    }
    return details.start !== undefined && details.end !== undefined
        ? `${details.referenceSequence}:${details.start}-${details.end}`
        : details.referenceSequence;
};

// "CD74::ROS1" (HGNC fusion notation)
export const getFusionName = (details: StructuralDetails) => details.fusionGenes?.join('::');

/**
 * The Explorer's label: a heading (the change type, or the fusion) and what follows it in parentheses, e.g.
 * "Copy number variation" + ["NC_000006.11:117609462-117747018", "1 copy"]
 */
export const getStructuralLabelParts = (details: StructuralDetails) => {
    if (details.fusionGenes) {
        return { heading: `${getFusionName(details)} gene fusion`, qualifiers: [] as string[] };
    }
    const qualifiers = [getStructuralLocation(details), getCopiesText(details)].filter((part): part is string => Boolean(part));
    return { heading: getChangeTypeLabel(details), qualifiers };
};

// The Explorer's label as text, e.g. "Copy number variation (NC_000006.11:117609462-117747018) (1 copy)"
export const getStructuralDescription = (details: StructuralDetails) => {
    const { heading, qualifiers } = getStructuralLabelParts(details);
    return [heading, ...qualifiers.map((qualifier) => `(${qualifier})`)].join(' ');
};

// The report card's change text, beside the gene: "Copy number variation (1 copy)", "Gene fusion"
export const getStructuralChange = (details: StructuralDetails) => {
    if (details.fusionGenes) {
        return 'Gene fusion';
    }
    const copies = getCopiesText(details);
    return copies ? `${getChangeTypeLabel(details)} (${copies})` : getChangeTypeLabel(details);
};

// Short label for pathway diagrams, e.g. "L858R", "CNV (12 copies)", "CD74::ROS1", "InDel"
export const getVariantShortLabel = (variant: Variant) => {
    const proteinChange = variant.molecularConsequences?.[0]?.proteinChange;

    if (proteinChange) {
        const change = proteinChange.split(':').pop()?.replace(/^p\./, '').replace(/[()]/g, '') ?? proteinChange;
        return change.replace(/[A-Z][a-z]{2}/g, (code) => AMINO_ACIDS[code] ?? code);
    }

    if (variant.structural) {
        if (variant.structural.fusionGenes) {
            return getFusionName(variant.structural) ?? 'Fusion';
        }
        const changeType = /copy_?number/i.test(variant.structural.changeType ?? '') ? 'CNV' : getChangeTypeLabel(variant.structural);
        const copies = getCopiesText(variant.structural);
        return copies ? `${changeType} (${copies})` : changeType;
    }

    return getSimpleVariantLabel(variant.variant);
};
