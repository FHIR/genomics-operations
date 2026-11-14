import { StructuralDetails, Variant } from '@/types/variants';
import { getStructuralDescription } from '@/lib/variantDisplay';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import { NormalizedVariant } from '@/services/reportedVariantLookupService';
import { ReportedVariantRecord } from '@/lib/reportedVariants';

// Layering the patient's reported variants onto the Genomics Explorer's results:
// the union of the server's variants and the reported ones,
// within the searched genes and ranges, with the reported ones flagged and their implications credited.
// Structural variants match and place by their recorded fields; fusions (not on the server) place by their partner genes.

export interface SearchedTerm {
  term: string;
  // The genomic range the term covers; unknown until resolved
  range?: string;
}

interface Spdi {
  accession: string;
  position: number; // 0-based
  deleted: string;
}

interface GenomicRange {
  accession: string;
  start: number;
  end: number;
}

function parseSpdi(value: string | undefined): Spdi | undefined {
  const match = value?.match(/^(NC_\d+\.\d+):(\d+):([A-Za-z]*):([A-Za-z]*)$/);
  return match ? { accession: match[1], position: Number(match[2]), deleted: match[3] } : undefined;
}

function parseRange(value: string | undefined): GenomicRange | undefined {
  const match = value?.match(/^(NC_\d+\.\d+):(\d+)-(\d+)$/);
  return match ? { accession: match[1], start: Number(match[2]), end: Number(match[3]) } : undefined;
}

// NC_000007.14 and NC_000007.13 are chromosome 7 in builds 38 and 37
function chromosome(accession: string) {
  return accession.split('.')[0];
}

// A simple substitution in gHGVS is an SPDI by arithmetic: g.POS REF>ALT is POS-1:REF:ALT
function substitutionSpdi(gHgvs: string | undefined) {
  const match = gHgvs?.match(/^(NC_\d+\.\d+):g\.(\d+)([ACGTN])>([ACGTN])$/);
  return match ? `${match[1]}:${Number(match[2]) - 1}:${match[3]}:${match[4]}` : undefined;
}

// Every SPDI known for a reported variant: its own, from its gHGVS, and from normalize-variant once looked up
function candidateSpdis(record: ReportedVariantRecord, normalized?: NormalizedVariant | null) {
  return [record.spdi, substitutionSpdi(record.gHgvs), normalized?.b37SPDI, normalized?.b38SPDI]
    .filter((spdi): spdi is string => Boolean(spdi));
}

// Same change type and location, and the same copy number when both give one (exact: overlap could pair different segments)
function sameStructuralVariant(reported: StructuralDetails, server: StructuralDetails) {
  return !reported.fusionGenes && !server.fusionGenes
    && reported.changeType === server.changeType
    && reported.referenceSequence === server.referenceSequence
    && reported.start !== undefined && reported.start === server.start
    && reported.end !== undefined && reported.end === server.end
    && (reported.copyNumber === undefined || server.copyNumber === undefined || reported.copyNumber === server.copyNumber);
}

/**
 * Whether a server row is this reported variant: by SPDI, by c.HGVS (same transcript and version), or by a
 * normalized SPDI; a structural variant by its recorded fields. Never by protein text.
 */
export function matchesRow(record: ReportedVariantRecord, row: Variant, normalized?: NormalizedVariant | null) {
  if (row.isReportOnly) {
    return false;
  }
  if (record.structural || row.variantType === 'structural') {
    return Boolean(record.structural && row.structural && sameStructuralVariant(record.structural, row.structural));
  }
  if (candidateSpdis(record, normalized).includes(row.variant)) {
    return true;
  }
  return Boolean(record.cHgvs) && (row.molecularConsequences ?? []).some((consequence) => consequence.dnaChange === record.cHgvs);
}

// Whether two spans on the same reference sequence share any position
function overlaps(start: number, end: number, target: GenomicRange) {
  return start <= target.end && end >= target.start;
}

// A structural variant (not a fusion) overlaps the range: true or false in the same build, false on another
// chromosome, undefined in another build of the same chromosome (no lookup can tell)
function structuralWithinRange(structural: StructuralDetails, target: GenomicRange) {
  const { referenceSequence, start, end } = structural;
  if (!referenceSequence || start === undefined || end === undefined) {
    return undefined;
  }
  if (referenceSequence === target.accession) {
    return overlaps(start, end, target);
  }
  return chromosome(referenceSequence) === chromosome(target.accession) ? undefined : false;
}

/**
 * Whether a fusion belongs under a searched term: the term is one of its partner genes, or a range that one
 * partner's coordinates (either build) overlap. False until those coordinates are known.
 */
export function fusionInTerm(fusionGenes: string[], term: SearchedTerm, geneRanges: Map<string, string[]>) {
  if (fusionGenes.includes(term.term)) {
    return true;
  }
  const target = parseRange(term.term);
  return Boolean(target) && fusionGenes.some((gene) => (geneRanges.get(gene) ?? [])
    .map(parseRange)
    .some((partner) => partner !== undefined && target !== undefined && partner.accession === target.accession
      && overlaps(partner.start, partner.end, target)));
}

/**
 * Whether the variant's genomic span lies within a range: true or false when a known SPDI is in the range's
 * build, false when it's on another chromosome, undefined when only a lookup can tell. A structural variant
 * needs only to overlap the range; a fusion is placed by fusionInTerm instead.
 */
export function withinRange(record: ReportedVariantRecord, range: string | undefined, normalized?: NormalizedVariant | null) {
  const target = parseRange(range);
  if (!target) {
    return undefined;
  }
  if (record.structural) {
    return record.structural.fusionGenes ? false : structuralWithinRange(record.structural, target);
  }
  const spdis = candidateSpdis(record, normalized).map(parseSpdi).filter((spdi): spdi is Spdi => spdi !== undefined);
  const sameBuild = spdis.find((spdi) => spdi.accession === target.accession);
  if (sameBuild) {
    return sameBuild.position >= target.start && sameBuild.position + sameBuild.deleted.length <= target.end;
  }
  if (spdis.length > 0 && spdis.every((spdi) => chromosome(spdi.accession) !== chromosome(target.accession))) {
    return false;
  }
  return undefined;
}

/**
 * A reported variant needs a normalize-variant lookup for a term when the cheap comparisons can't settle it: no row
 * matches, and its location in the range's build is unknown.
 */
export function needsNormalization(record: ReportedVariantRecord, term: SearchedTerm, rows: Variant[]) {
  // normalize-variant takes SPDI or HGVS, which structural variants and fusions don't have
  return !record.structural && !rows.some((row) => matchesRow(record, row)) && term.range !== undefined && withinRange(record, term.range) === undefined;
}

function sameImplication(left: ProcessedTxImplication, right: ProcessedTxImplication) {
  const identity = (implication: ProcessedTxImplication) => [
    implication.medication, implication.therapeuticImplication, implication.phenotypicContext, implication.evidenceLevel,
  ].map((part) => (part ?? '').trim().toLowerCase()).join('|');
  return identity(left) === identity(right);
}

function reportImplications(record: ReportedVariantRecord): ProcessedTxImplication[] {
  return record.implications.map((implication) => ({
    phenotypicContext: implication.phenotype ?? '',
    phenotypeCodings: implication.phenotypeCodings?.flatMap(({ system, code, display }) => (code ? [{ system, code, display }] : [])),
    evidenceLevel: implication.evidence ?? '',
    medication: implication.drug ?? '',
    therapeuticImplication: implication.implication ?? '',
    hyperlink: '',
    reportSources: implication.reports,
  }));
}

// The report's implications first, each credited to its reports; one the server also has keeps the server's
// CIViC link. Then the server's other implications, as before.
function mergeImplications(server: ProcessedTxImplication[], fromReports: ProcessedTxImplication[]) {
  const used = new Set<number>();
  const merged = fromReports.map((reported) => {
    const index = server.findIndex((candidate, position) => !used.has(position) && sameImplication(candidate, reported));
    if (index < 0) {
      return reported;
    }
    used.add(index);
    return { ...server[index], reportSources: reported.reportSources };
  });
  return [...merged, ...server.filter((_, position) => !used.has(position))];
}

function uniqueCredits(records: ReportedVariantRecord[]) {
  const seen = new Set<string>();
  return records.flatMap((record) => record.reports).filter((credit) => {
    const id = `${credit.title}|${credit.date ?? ''}`;
    return !seen.has(id) && Boolean(seen.add(id));
  });
}

function reportOnlyRow(record: ReportedVariantRecord, term: SearchedTerm, normalized?: NormalizedVariant | null): Variant {
  const shared = {
    id: `reported|${term.term}|${record.key}`,
    range: term.term,
    resolvedRange: term.range,
    genomicSourceClass: record.origin,
    variantAlleleFrequency: record.vaf,
    dxImplications: [],
    txImplications: reportImplications(record),
    isLoading: false,
    reportedIn: record.reports,
    isReportOnly: true,
  };
  if (record.structural) {
    return { ...shared, variant: getStructuralDescription(record.structural), variantType: 'structural', structural: record.structural, molecularConsequences: [] };
  }
  const range = parseRange(term.range);
  const spdis = candidateSpdis(record, normalized);
  // Shown in the range's build when known, like the server's rows
  const variant = spdis.find((spdi) => parseSpdi(spdi)?.accession === range?.accession) ?? spdis[0] ?? record.gHgvs ?? record.cHgvs ?? record.key;
  const consequences: MolecularConsequence[] = record.pHgvs || record.cHgvs
    ? [{ featureConsequence: '', interpretation: '', impact: '', hyperlink: '', proteinChange: record.pHgvs, dnaChange: record.cHgvs }]
    : [];
  return { ...shared, variant, variantType: 'simple', molecularConsequences: consequences };
}

/**
 * The union: each searched term's server rows, reported ones flagged with their reports and implications
 * merged, then the reported variants within the term's range that the server doesn't have, and the reported
 * fusions with a partner in the term (geneRanges: each fusion partner's coordinates, for range terms).
 */
export function layerReportedVariants(
  serverRows: Variant[],
  terms: SearchedTerm[],
  records: ReportedVariantRecord[],
  normalized: Map<string, NormalizedVariant | null>,
  geneRanges: Map<string, string[]>,
): Variant[] {
  if (records.length === 0) {
    return serverRows;
  }
  const output: Variant[] = [];
  const termNames = new Set(terms.map((term) => term.term));
  for (const term of terms) {
    const rows = serverRows.filter((row) => row.range === term.term);
    const matchedKeys = new Set<string>();
    for (const row of rows) {
      const matches = records.filter((record) => matchesRow(record, row, normalized.get(record.key)));
      matches.forEach((record) => matchedKeys.add(record.key));
      output.push(matches.length === 0 ? row : {
        ...row,
        reportedIn: uniqueCredits(matches),
        txImplications: mergeImplications(row.txImplications ?? [], matches.flatMap(reportImplications)),
      });
    }
    for (const record of records) {
      const fusionGenes = record.structural?.fusionGenes;
      const inTerm = fusionGenes
        ? fusionInTerm(fusionGenes, term, geneRanges)
        : withinRange(record, term.range, normalized.get(record.key)) === true;
      if (!matchedKeys.has(record.key) && inTerm) {
        output.push(reportOnlyRow(record, term, normalized.get(record.key)));
      }
    }
  }
  // Rows from a term no longer listed are kept as they were
  return [...output, ...serverRows.filter((row) => !termNames.has(row.range))];
}
