import { Variant } from '@/types/variants';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import { NormalizedVariant } from '@/services/reportedVariantLookupService';
import { ReportedVariantRecord } from '@/lib/reportedVariants';

// Layering the patient's reported variants onto the Genomics Explorer's results:
// the union of the server's variants and the reported ones,
// within the searched genes and ranges, with the reported ones flagged and their implications credited

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

/**
 * Whether a server row is this reported variant: by SPDI, by c.HGVS (same transcript and version), or by a
 * normalized SPDI. Never by protein text; structural-variant rows never match.
 */
export function matchesRow(record: ReportedVariantRecord, row: Variant, normalized?: NormalizedVariant | null) {
  if (row.variantType === 'structural' || row.isReportOnly) {
    return false;
  }
  if (candidateSpdis(record, normalized).includes(row.variant)) {
    return true;
  }
  return Boolean(record.cHgvs) && (row.molecularConsequences ?? []).some((consequence) => consequence.dnaChange === record.cHgvs);
}

/**
 * Whether the variant's genomic span lies within a range: true or false when a known SPDI is in the range's
 * build, false when it's on another chromosome, undefined when only a lookup can tell.
 */
export function withinRange(record: ReportedVariantRecord, range: string | undefined, normalized?: NormalizedVariant | null) {
  const target = parseRange(range);
  if (!target) {
    return undefined;
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
  return !rows.some((row) => matchesRow(record, row)) && term.range !== undefined && withinRange(record, term.range) === undefined;
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
  const range = parseRange(term.range);
  const spdis = candidateSpdis(record, normalized);
  // Shown in the range's build when known, like the server's rows
  const variant = spdis.find((spdi) => parseSpdi(spdi)?.accession === range?.accession) ?? spdis[0] ?? record.gHgvs ?? record.cHgvs ?? record.key;
  const consequences: MolecularConsequence[] = record.pHgvs || record.cHgvs
    ? [{ featureConsequence: '', interpretation: '', impact: '', hyperlink: '', proteinChange: record.pHgvs, dnaChange: record.cHgvs }]
    : [];
  return {
    id: `reported|${term.term}|${record.key}`,
    range: term.term,
    resolvedRange: term.range,
    variant,
    variantType: 'simple',
    genomicSourceClass: record.origin,
    variantAlleleFrequency: record.vaf,
    dxImplications: [],
    txImplications: reportImplications(record),
    molecularConsequences: consequences,
    isLoading: false,
    reportedIn: record.reports,
    isReportOnly: true,
  };
}

/**
 * The union: each searched term's server rows, reported ones flagged with their reports and implications
 * merged, then the reported variants within the term's range that the server doesn't have.
 */
export function layerReportedVariants(
  serverRows: Variant[],
  terms: SearchedTerm[],
  records: ReportedVariantRecord[],
  normalized: Map<string, NormalizedVariant | null>,
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
      if (!matchedKeys.has(record.key) && withinRange(record, term.range, normalized.get(record.key)) === true) {
        output.push(reportOnlyRow(record, term, normalized.get(record.key)));
      }
    }
  }
  // Rows from a term no longer listed are kept as they were
  return [...output, ...serverRows.filter((row) => !termNames.has(row.range))];
}
