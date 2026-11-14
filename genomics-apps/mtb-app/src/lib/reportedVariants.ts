import { Coding, DiagnosticReport, FhirResource, Observation, Reference } from '@/types/fhir';
import { ReportCredit, StructuralDetails } from '@/types/variants';
import { getFusionName, getStructuralChange, getStructuralLocation } from '@/lib/variantDisplay';
import { conceptLabel, formatDate, LOINC } from '@/lib/clinicalSummary';
import { TxComponentCodes } from '@/utils/ComponentCodes';

// The Reported variants card: the variants and therapeutic
// implications that came into the record from the lab's genetic reports, grouped by report

const GR_TBD = 'http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs';
const CLINVAR = 'http://www.ncbi.nlm.nih.gov/clinvar';
const GENETIC_ANALYSIS_REPORT = '51969-4';
const VARIANT_CODE = '69548-6'; // Genetic variant assessment
const REGION_STUDIED_CODE = '53041-0'; // DNA region of interest panel
const PRESENT = 'LA9633-4';
const COMPONENT = {
  gene: '48018-6',
  cHgvs: '48004-6',
  pHgvs: '48005-3',
  discreteVariant: '81252-9',
  gHgvs: '81290-9',
  origin: '48002-0',
  changeType: '48019-4',
  copyNumber: '82155-3',
  referenceSequence: '48013-7',
  outerRange: '81301-4',
  innerRange: '81302-2',
  altAllele: '69551-0',
  molecularConsequence: 'molecular-consequence',
};
const GENE_FUSION = 'SO:0001565';
const EVIDENCE_ORDER = ['A', 'B', 'C', 'D', 'E'];

export type EvidenceLevel = 'A' | 'B' | 'C' | 'D' | 'E';

export interface ReportedImplication {
  key: string;
  // Either may be absent when not recorded; nothing is shown in its place
  drug?: string;
  implication?: string;
  // The level as recorded, e.g. "B (Clinical evidence)", and its letter when it starts with one
  evidence?: string;
  level?: EvidenceLevel;
  phenotype?: string;
  phenotypeCodings?: Coding[]; // The cancer types the evidence applies to, for matching the patient's tumor type
}

export interface ReportedVariant {
  key: string;
  // Absent when the report doesn't name the gene; a fusion's partners joined, e.g. "CD74::ROS1"
  gene?: string;
  // The change without its transcript prefix, e.g. "c.3549C>A (p.Tyr1183Ter)"; absent when not recorded
  change?: string;
  // Full HGVS, genomic location and ClinVar id, for the hover
  details: string[];
  // Genomic source class, e.g. "Germline"
  origin?: string;
  implications: ReportedImplication[];
  // Identity, for matching against the genomics server's variants (by SPDI, then c.HGVS, then normalized SPDI)
  spdi?: string; // NCBI variation SPDI, when the report has one
  cHgvs?: string;
  pHgvs?: string;
  gHgvs?: string;
  vaf?: number;
  // A structural variant's or fusion's fields (no HGVS); absent for a simple variant
  structural?: StructuralDetails;
}

export interface GeneticReportGroup {
  key: string;
  title: string;
  date?: string;
  // Specimen and performer, for the header's hover
  details: string[];
  variants: ReportedVariant[];
  // Shown when the report has no structured variants
  conclusion?: string;
}

export interface ReportedVariantsData {
  groups: GeneticReportGroup[];
  // Implications whose variant isn't a reported variant in this record
  otherImplications: ReportedImplication[];
  // Every variant, region-studied and implication Observation, so other sections leave them out
  genomicObservations: Observation[];
}

export type { ReportCredit };

// One reported variant across all the patient's reports: a variant in two reports (L2345's EGFR L858R) is one
// record crediting both, and an implication in both counts once
export interface ReportedVariantRecord extends Omit<ReportedVariant, 'implications'> {
  reports: ReportCredit[];
  implications: (ReportedImplication & { reports: ReportCredit[] })[];
}

type Component = NonNullable<Observation['component']>[number];

function key(resource: FhirResource) {
  return `${resource.resourceType}/${resource.id}`;
}

function referenceKey(reference: Reference | undefined) {
  const parts = reference?.reference?.split('/') ?? [];
  return parts.length >= 2 ? parts.slice(-2).join('/') : undefined;
}

function codeIs(concept: { coding?: { system?: string; code?: string }[] } | undefined, system: string, code: string) {
  return concept?.coding?.some((coding) => coding.system === system && coding.code === code) ?? false;
}

export function isGeneticVariant(observation: Observation) {
  return codeIs(observation.code, LOINC, VARIANT_CODE);
}

function isRegionStudied(observation: Observation) {
  return codeIs(observation.code, LOINC, REGION_STUDIED_CODE);
}

export function isTherapeuticImplication(observation: Observation) {
  return codeIs(observation.code, GR_TBD, 'therapeutic-implication');
}

function isGeneticReport(report: DiagnosticReport) {
  return codeIs(report.code, LOINC, GENETIC_ANALYSIS_REPORT);
}

function component(observation: Observation, code: string): Component | undefined {
  return observation.component?.find((candidate) => candidate.code?.coding?.some((coding) => coding.code === code));
}

function componentCoding(observation: Observation, code: string) {
  return component(observation, code)?.valueCodeableConcept?.coding?.[0];
}

function componentText(observation: Observation, code: string) {
  const value = component(observation, code);
  return conceptLabel(value?.valueCodeableConcept) ?? value?.valueString;
}

// Every gene component's symbol, in recorded order (a fusion has one per partner)
function geneSymbols(observation: Observation) {
  return (observation.component ?? [])
    .filter((candidate) => candidate.code?.coding?.some((coding) => coding.code === COMPONENT.gene))
    .map((candidate) => conceptLabel(candidate.valueCodeableConcept) ?? candidate.valueString)
    .filter((gene): gene is string => Boolean(gene));
}

function isFusion(observation: Observation) {
  return (observation.component ?? []).some((candidate) => candidate.code?.coding?.some((coding) => coding.code === COMPONENT.molecularConsequence)
    && candidate.valueCodeableConcept?.coding?.some((coding) => coding.code === GENE_FUSION));
}

// A structural variant's fields (DNA change type, location, copy number) or a fusion's partner genes;
// undefined for a simple variant
function structuralDetails(observation: Observation): StructuralDetails | undefined {
  if (isFusion(observation)) {
    return { fusionGenes: geneSymbols(observation) };
  }
  const changeType = componentCoding(observation, COMPONENT.changeType)?.display;
  if (!changeType) {
    return undefined;
  }
  const outer = component(observation, COMPONENT.outerRange)?.valueRange;
  const inner = component(observation, COMPONENT.innerRange)?.valueRange;
  const range = outer?.low?.value !== undefined && outer?.high?.value !== undefined ? outer : inner;
  const insertedSequence = component(observation, COMPONENT.altAllele)?.valueString;
  return {
    changeType,
    referenceSequence: componentCoding(observation, COMPONENT.referenceSequence)?.code,
    start: range?.low?.value,
    end: range?.high?.value,
    copyNumber: component(observation, COMPONENT.copyNumber)?.valueQuantity?.value,
    insertedSequence: insertedSequence && /^[ACGTN]+$/i.test(insertedSequence) ? insertedSequence : undefined,
  };
}

// "NM_006015.6:c.5161C>T" becomes "c.5161C>T"; "NP_006006.3:p.(Val2169CysfsTer56)" becomes "p.(Val2169CysfsTer56)"
function withoutAccession(hgvs: string | undefined) {
  return hgvs?.replace(/^[^:]*:/, '');
}

function variantChange(observation: Observation) {
  const cHgvs = withoutAccession(componentCoding(observation, COMPONENT.cHgvs)?.code);
  const pHgvs = withoutAccession(componentCoding(observation, COMPONENT.pHgvs)?.code);
  if (cHgvs) {
    return pHgvs ? `${cHgvs} (${pHgvs})` : cHgvs;
  }
  // A ClinVar name such as "NM_024675.3(PALB2):c.3549C>A (p.Tyr1183Ter)" loses its transcript and gene
  const discrete = componentCoding(observation, COMPONENT.discreteVariant);
  if (discrete?.system === CLINVAR && discrete.display) {
    return discrete.display.replace(/^[^:]*\):/, '');
  }
  return componentCoding(observation, COMPONENT.gHgvs)?.code ?? discrete?.code;
}

function variantDetails(observation: Observation, structural: StructuralDetails | undefined) {
  const discrete = componentCoding(observation, COMPONENT.discreteVariant);
  const notes = ((observation.note as { text?: string }[] | undefined) ?? []).map((note) => note.text);
  return [
    structural && getStructuralLocation(structural),
    structural?.insertedSequence && `Inserted sequence: ${structural.insertedSequence}`,
    // A fusion's note, e.g. "CD74-ROS1 fusion; in-frame"
    ...(structural?.fusionGenes ? notes : []),
    componentCoding(observation, COMPONENT.cHgvs)?.code,
    componentCoding(observation, COMPONENT.pHgvs)?.code,
    discrete?.system === CLINVAR ? discrete.display : undefined,
    componentCoding(observation, COMPONENT.gHgvs)?.code,
    discrete && discrete.system !== CLINVAR ? discrete.code : undefined,
    discrete?.system === CLINVAR ? `ClinVar ${discrete.code}` : undefined,
  ].filter((detail): detail is string => Boolean(detail));
}

// A ClinVar name such as "NM_024675.3(PALB2):c.3549C>A (p.Tyr1183Ter)" gives the c.HGVS ("NM_024675.3:c.3549C>A")
// and p.HGVS ("p.Tyr1183Ter") when the report has no c./p. components (MRN1234's PALB2)
function clinVarHgvs(observation: Observation) {
  const discrete = componentCoding(observation, COMPONENT.discreteVariant);
  const match = discrete?.system === CLINVAR ? discrete.display?.match(/^([^(\s]+)\([^)]*\):(c\.[^\s]+)(?:\s+\((p\.[^)]+)\))?$/) : undefined;
  return match ? { cHgvs: `${match[1]}:${match[2]}`, pHgvs: match[3] } : undefined;
}

// The Discrete genetic variant component when it's an NCBI variation SPDI (not a ClinVar id)
function discreteSpdi(observation: Observation) {
  const discrete = componentCoding(observation, COMPONENT.discreteVariant);
  return discrete?.system?.includes('ncbi.nlm.nih.gov/variation') ? discrete.code : undefined;
}

function capitalized(text: string | undefined) {
  return text ? text.charAt(0).toUpperCase() + text.slice(1) : undefined;
}

function evidenceRank(level: EvidenceLevel | undefined) {
  return level ? EVIDENCE_ORDER.indexOf(level) : EVIDENCE_ORDER.length;
}

// A combination therapy has one medication-assessed component per drug: "Docetaxel + Selumetinib"
function drugs(observation: Observation) {
  const names = (observation.component ?? [])
    .filter((candidate) => candidate.code?.coding?.some((coding) => coding.code === TxComponentCodes.medication))
    .map((candidate) => capitalized(conceptLabel(candidate.valueCodeableConcept) ?? candidate.valueString))
    .filter((name): name is string => Boolean(name));
  return names.length > 0 ? names.join(' + ') : undefined;
}

// An implication can apply to several cancer types, one component each
function phenotypeComponents(observation: Observation) {
  return (observation.component ?? [])
    .filter((candidate) => candidate.code?.coding?.some((coding) => coding.code === TxComponentCodes.phenotypicContext));
}

function phenotypes(observation: Observation) {
  const names = phenotypeComponents(observation)
    .map((candidate) => conceptLabel(candidate.valueCodeableConcept) ?? candidate.valueString)
    .filter((name): name is string => Boolean(name));
  return names.length > 0 ? names.join('; ') : undefined;
}

function toImplication(observation: Observation): ReportedImplication {
  const evidence = componentText(observation, TxComponentCodes.evidenceLevel);
  const letter = evidence?.match(/^([A-E])\b/)?.[1] as EvidenceLevel | undefined;
  return {
    key: key(observation),
    drug: drugs(observation),
    implication: componentText(observation, TxComponentCodes.therapeuticImplication),
    evidence,
    level: letter,
    phenotype: phenotypes(observation),
    phenotypeCodings: phenotypeComponents(observation).flatMap((candidate) => candidate.valueCodeableConcept?.coding ?? []),
  };
}

const byEvidenceThenDrug = (left: ReportedImplication, right: ReportedImplication) =>
  evidenceRank(left.level) - evidenceRank(right.level) || (left.drug ?? '').localeCompare(right.drug ?? '');

// Best implication's evidence first (variants with none last), then gene, then change
function byBestEvidenceThenGene(left: ReportedVariant, right: ReportedVariant) {
  const best = (variant: ReportedVariant) => Math.min(...variant.implications.map((implication) => evidenceRank(implication.level)), EVIDENCE_ORDER.length + 1);
  return best(left) - best(right) || (left.gene ?? '').localeCompare(right.gene ?? '') || (left.change ?? '').localeCompare(right.change ?? '');
}

function reportDetails(report: DiagnosticReport, index: Map<string, FhirResource>) {
  const specimens = ((report.specimen as Reference[] | undefined) ?? []).map((reference) => {
    const specimen = index.get(referenceKey(reference) ?? '');
    const note = (specimen?.note as { text?: string }[] | undefined)?.[0]?.text;
    return reference.display ?? note ?? conceptLabel(specimen?.type as Parameters<typeof conceptLabel>[0]);
  });
  const performers = ((report.performer as Reference[] | undefined) ?? []).map((reference) => reference.display);
  return [...specimens, ...performers].filter((detail): detail is string => Boolean(detail));
}

export function buildReportedVariants(resources: FhirResource[], index: Map<string, FhirResource>): ReportedVariantsData {
  const observations = resources.filter((resource): resource is Observation => resource.resourceType === 'Observation');
  const variantObservations = observations.filter((observation) => isGeneticVariant(observation)
    && (!observation.valueCodeableConcept || codeIs(observation.valueCodeableConcept, LOINC, PRESENT)));
  const implicationObservations = observations.filter(isTherapeuticImplication);
  const genomicObservations = observations.filter((observation) => isGeneticVariant(observation) || isRegionStudied(observation) || isTherapeuticImplication(observation));

  // Each variant with the implications derived from it
  const variants = new Map<string, ReportedVariant>(variantObservations.map((observation) => {
    const structural = structuralDetails(observation);
    return [key(observation), {
      key: key(observation),
      gene: (structural?.fusionGenes && getFusionName(structural)) || componentText(observation, COMPONENT.gene),
      change: variantChange(observation) ?? (structural && getStructuralChange(structural)),
      details: variantDetails(observation, structural),
      origin: capitalized(componentText(observation, COMPONENT.origin)),
      implications: [],
      spdi: discreteSpdi(observation),
      cHgvs: componentCoding(observation, COMPONENT.cHgvs)?.code ?? clinVarHgvs(observation)?.cHgvs,
      pHgvs: componentCoding(observation, COMPONENT.pHgvs)?.code ?? clinVarHgvs(observation)?.pHgvs,
      gHgvs: componentCoding(observation, COMPONENT.gHgvs)?.code,
      vaf: component(observation, '81258-6')?.valueQuantity?.value,
      structural,
    }];
  }));
  const otherImplications: ReportedImplication[] = [];
  for (const observation of implicationObservations) {
    const implication = toImplication(observation);
    const targets = ((observation.derivedFrom as Reference[] | undefined) ?? [])
      .map((reference) => variants.get(referenceKey(reference) ?? ''))
      .filter((variant): variant is ReportedVariant => variant !== undefined);
    if (targets.length === 0) {
      otherImplications.push(implication);
    }
    targets.forEach((variant) => variant.implications.push(implication));
  }
  variants.forEach((variant) => variant.implications.sort(byEvidenceThenDrug));

  // One group per genetic report, newest first; variants no report lists go in a final group
  const reports = resources.filter((resource): resource is DiagnosticReport => resource.resourceType === 'DiagnosticReport')
    .filter((report) => isGeneticReport(report) || (report.result ?? []).some((reference) => variants.has(referenceKey(reference) ?? '')));
  const dateOf = (report: DiagnosticReport) => report.effectiveDateTime ?? report.effectivePeriod?.start ?? report.issued;
  const listed = new Set<string>();
  const groups: GeneticReportGroup[] = reports
    .sort((left, right) => (dateOf(right) ?? '').localeCompare(dateOf(left) ?? ''))
    .map((report) => {
      const reportVariants = (report.result ?? [])
        .map((reference) => variants.get(referenceKey(reference) ?? ''))
        .filter((variant): variant is ReportedVariant => variant !== undefined);
      reportVariants.forEach((variant) => listed.add(variant.key));
      return {
        key: key(report),
        title: conceptLabel(report.code) ?? 'Genetic report',
        date: formatDate(dateOf(report)),
        details: reportDetails(report, index),
        variants: reportVariants.sort(byBestEvidenceThenGene),
        conclusion: reportVariants.length === 0 ? report.conclusion : undefined,
      };
    });
  const unlisted = [...variants.values()].filter((variant) => !listed.has(variant.key));
  if (unlisted.length > 0) {
    groups.push({ key: 'unlisted', title: 'Variants not in a report', details: [], variants: unlisted.sort(byBestEvidenceThenGene) });
  }
  otherImplications.sort(byEvidenceThenDrug);
  return { groups, otherImplications, genomicObservations };
}

function implicationIdentity(implication: ReportedImplication) {
  return [implication.drug, implication.implication, implication.phenotype, implication.evidence]
    .map((part) => (part ?? '').toLowerCase()).join('|');
}

// The same structural variant in two reports is one record: same gene, change type, location and copy number
// (genes sharing one copy-number segment stay separate); a fusion, by its partner genes
function structuralIdentity(structural: StructuralDetails, gene?: string) {
  return structural.fusionGenes
    ? `fusion|${structural.fusionGenes.join('::')}`
    : ['structural', gene, structural.changeType, structural.referenceSequence, structural.start, structural.end, structural.copyNumber].join('|');
}

/**
 * The patient's reported variants, one record per distinct variant (by SPDI, else gHGVS, else c.HGVS; a
 * structural variant or fusion by its recorded fields), each
 * crediting the reports it was in; implications that more than one report gives count once.
 */
export function collectReportedVariants(data: ReportedVariantsData): ReportedVariantRecord[] {
  const records = new Map<string, ReportedVariantRecord>();
  for (const group of data.groups) {
    const credit: ReportCredit = { title: group.title, date: group.date };
    for (const variant of group.variants) {
      const identity = variant.spdi ?? variant.gHgvs ?? variant.cHgvs ?? (variant.structural && structuralIdentity(variant.structural, variant.gene)) ?? variant.key;
      let record = records.get(identity);
      if (!record) {
        record = { ...variant, reports: [], implications: [] };
        records.set(identity, record);
      }
      record.reports.push(credit);
      for (const implication of variant.implications) {
        const existing = record.implications.find((candidate) => implicationIdentity(candidate) === implicationIdentity(implication));
        if (existing) {
          existing.reports.push(credit);
        } else {
          record.implications.push({ ...implication, reports: [credit] });
        }
      }
    }
  }
  return [...records.values()];
}
