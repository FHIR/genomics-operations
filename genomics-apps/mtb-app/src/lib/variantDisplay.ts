import { Variant } from '@/types/variants';

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

// Short label for pathway diagrams, e.g. "L858R", "CNV (12 copies)", "InDel"
export const getVariantShortLabel = (variant: Variant) => {
    const proteinChange = variant.molecularConsequences?.[0]?.proteinChange;

    if (proteinChange) {
        const change = proteinChange.split(':').pop()?.replace(/^p\./, '').replace(/[()]/g, '') ?? proteinChange;
        return change.replace(/[A-Z][a-z]{2}/g, (code) => AMINO_ACIDS[code] ?? code);
    }

    if (variant.variantType === 'structural') {
        const rawChangeType = variant.variant.split(' ')[0];
        const changeType = /copy_?number/i.test(rawChangeType) ? 'CNV' : rawChangeType.replace(/_/g, ' ');
        const copies = /Copies:\s*([\d.]+)/.exec(variant.variant)?.[1];
        return copies ? `${changeType} (${copies} copies)` : changeType;
    }

    return getSimpleVariantLabel(variant.variant);
};
