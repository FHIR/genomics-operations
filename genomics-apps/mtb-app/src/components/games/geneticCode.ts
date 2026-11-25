// The standard genetic code, for the Mutation Hunter game

const BASES = ['T', 'C', 'A', 'G'] as const;
// Amino acids for codons in TTT, TTC, TTA, TTG, TCT … GGG order ('*' is a stop codon)
const AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';

export const CODON_TABLE: Record<string, string> = {};
BASES.forEach((first, i) => BASES.forEach((second, j) => BASES.forEach((third, k) => {
  CODON_TABLE[first + second + third] = AMINO_ACIDS[i * 16 + j * 4 + k];
})));

export const THREE_LETTER: Record<string, string> = {
  A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile',
  L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter',
};

export const DNA_BASES = ['A', 'C', 'G', 'T'];

export function translate(sequence: string): string[] {
  const protein: string[] = [];
  for (let i = 0; i + 3 <= sequence.length; i += 3) {
    protein.push(CODON_TABLE[sequence.slice(i, i + 3)]);
  }
  return protein;
}

export const SENSE_CODONS = Object.keys(CODON_TABLE).filter((codon) => CODON_TABLE[codon] !== '*' && codon !== 'ATG');
