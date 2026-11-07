# How a case synthesis is made

This folder holds the case syntheses shown on the mtb-app's **AI Synthesis** tab, one Markdown file per patient, named `<MRN>.md`. This README describes how they are written, so that they can be revised, and written again when a patient's data changes, the same way each time.

## 1. Purpose and scope

A synthesis brings one patient's clinical record and genomic findings together into a short brief for a molecular tumor board. It covers who the patient is and where they are in treatment, what the tumor's genomic findings mean, what's missing, and which treatment options and clinical trials the evidence supports.

- It is a **research demonstration** on synthetic patients and is **not for clinical use**. Every synthesis says so at the top.
- It reports what labels, knowledge bases and guidelines support for this patient, and may recommend for or against an option when the evidence supports doing so.
- For the demo, syntheses are written in advance, not generated live. The tab's Generate button shows a short progress display (`SynthesisView.tsx`) and then loads the stored file (`synthesisService.ts`). A live AI call could later replace `fetchSynthesis` without changing the tab.
- A patient with no file here shows "Synthesis not available for this patient in the demo."

Current files: `HCC1143.md`, `L2345.md`, `MRN1234.md`, `TCGA-DD-A1EH.md`.

## 2. Inputs

Everything in a synthesis should come from the following sources, the same data the app shows on its other tabs. Use only what would be available for a real patient: the scenario notes in `clinical-data/<MRN>/` describe how the demo data was built and are not an input.

| Input | Where it comes from | Used for |
|---|---|---|
| Clinical record | The patient's FHIR bundle, `clinical-data/<MRN>/<MRN>-transaction-bundle.json`, loaded to the clinical FHIR server (the Case Summary tab) | Case summary, current status, prior treatment, information gaps |
| Genetic reports | The bundle's genomic reports, with their reported variants and implications (the Reported variants card) | Genomic findings as reported by the lab |
| Genomic data | The genomics server's operations (the Genomics Explorer tab), for the genes in §2.1 | Genomic findings and their evidence |
| Co-occurrence evidence | `public/data/CIViC_Cooccurrence_KB.csv` | Combination therapy for co-occurring alterations |
| Cancer genes | `public/data/OncoKB_CancerGeneList.tsv` | Which genes to search for cancer drivers, and each gene's role (§2.5) |

### 2.1 Which genes to search

The presets named below are the Explorer's; the same genes can be searched by calling the operations in §2.4 directly.

**Minimum, for every synthesis:**

- **Reported genes:** the genes of the patient's reported variants (the Explorer's Reported genes preset).
- **Gene lists:** every gene list in `public/data/MTB_KB_GeneLists.csv` for the patient's diagnosis (the Explorer's Actionable Gene List and Extended Gene List presets). Hepatocellular carcinoma has only an Extended Gene List. The ACMG Secondary Findings list is not part of the minimum.
- **Always-search genes:** the tumor-agnostic list in §2.3, searched for every patient whatever the tumor type.

### 2.2 Searching beyond the minimum

The minimum is a floor, not a ceiling. Whoever writes the synthesis (person or AI) is empowered to query any part of the patient's genome, as often and as widely as they judge necessary, guided by the literature, labels, trials, guidelines and clinical reasoning. The operations in §2.4 can be called directly, for any range in the patient's genome build; no further permission is needed. There is no limit on the number of queries or how much of the genome they cover; the only stopping rule is the loop below. Every query goes in the search log (below), which records why it was made but is not a gate on making it.

**The loop.** Gathering findings and looking up evidence (§4, steps 2–8) repeat until they stop raising new questions:

1. Search the minimum genes.
2. Read the evidence for what was found and for the patient's disease setting: drug labels, CIViC, guidelines, the literature, and trials.
3. Note every genomic question the evidence raises. For example, a resistance mechanism the literature describes after the patient's current therapy, or a trial's genomic eligibility criteria.
4. Query the genome to answer each question. Add any new finding to the genomic findings and go back to step 2.
5. Stop when the evidence raises no question that hasn't been queried.

**What can be queried.** Not only whole genes:

| Question | How |
|---|---|
| Any variant in a gene | `get-feature-coordinates` for the gene's range, then the simple and structural variant searches |
| A specific variant (e.g. EGFR C797S) | `normalize-variant` to get its canonical form, then `$find-subject-specific-variants` |
| A region (e.g. MET exon 14 splice sites, the TERT promoter) | The variant searches over that range |
| A structural event (fusion partner, copy-number change) | `$find-subject-structural-intersecting-variants` over each partner gene |
| A class of variant (e.g. frameshifts) | Search the gene, then filter by `$find-subject-molecular-consequences` |

Some questions can't be answered from variants: MSI, TMB, protein expression (IHC), or fusions confirmed by RNA testing. Unless the record has them, they go under Information gaps.

**Not found vs. not tested.** A search that returns nothing doesn't by itself mean the finding is absent. Call `$find-study-metadata` with the same `ranges` as the search. For each of the patient's tests it returns:

- the test date, specimen and **genome build** (search with coordinates in that build)
- the **variant types** (`dnaChangeType`) the test can detect. A test without CNV or structural variant types can't rule out an amplification or fusion.
- **`regionStudied`**: the parts of the range the test studied
- **`uncallableRegions`**: the parts of the range where the test couldn't make a call

A finding is **absent** only if its location is within a studied region, not within an uncallable region, and covered by a test that detects that variant type. Otherwise it's **not tested** (or not callable), and goes under Information gaps. For example, a panel that studied only MET's exons can rule out an exon 14 point mutation but not an intronic splice-site variant outside the studied regions.

**Search log.** Every lookup is recorded in a companion file, `<MRN>.searches.md`, next to the synthesis: genome queries, and also label, CIViC, ClinVar, trial and literature lookups. The app doesn't show it. Each entry gives:

- what was queried (gene, variant, or range and build; or the source and search terms)
- why (the minimum list, or the label, trial, guideline or paper that raised it)
- the result, and the date

The References section's "Genes searched" line summarizes the log.

### 2.3 Always-search genes (tumor-agnostic)

The loop in §2.2 starts from the patient's tumor type and findings, so it can miss an option that depends only on a genomic finding. For example, suppose a trial enrolls any tumor type with a CFTR frameshift. Nothing about a lung cancer patient leads to CFTR. Disease-based trial searches won't surface it, and CFTR isn't on any gene list. The patient's CFTR frameshift would only be found if CFTR were searched for everyone.

So these options are found from the other direction, and their genes are kept on a list searched for every patient:

1. **Find tumor-agnostic options:** FDA tissue-agnostic approvals, and recruiting ClinicalTrials.gov trials open to any solid tumor (or any cancer) whose eligibility rests on a genomic finding (basket and biomarker-driven trials).
2. **Pull out the genes and what's required:** gene, variant class (e.g. "frameshift", "fusion", "V600E"), and origin if specified (somatic or germline).
3. **Add them to the list below**, each with its source and the date checked. Remove entries whose approval is withdrawn or whose trial has stopped recruiting.
4. **When writing a synthesis**, search every gene on the list. A match counts only if it meets the entry's variant class and origin. A germline CFTR variant doesn't meet a trial for somatic frameshifts.

The list is refreshed before a batch of syntheses is written. A tumor-agnostic option found while working on one patient (in a trial search, say) is added to the list then, so it is searched for every later patient as well.

| Gene | Required finding | Source | Checked |
|---|---|---|---|
| NTRK1, NTRK2, NTRK3 | Fusion | FDA labels: Vitrakvi (larotrectinib) 2026-05-22, Rozlytrek (entrectinib) 2025-12-04, Augtyro (repotrectinib) 2024-06-13 | 2026-11-06 |
| BRAF | V600E | FDA labels: Tafinlar (dabrafenib) 2026-05-20 with Mekinist (trametinib) 2026-05-07 | 2026-11-06 |
| RET | Fusion | FDA label: Retevmo (selpercatinib) 2026-07-14 | 2026-11-06 |
| MLH1, MSH2, MSH6, PMS2 | Supports MSI-H/dMMR, which the labels require to be established by an FDA-authorized test | FDA labels: Keytruda (pembrolizumab, IV) 2026-07-10, Jemperli (dostarlimab) 2025-09-02 | 2026-11-06 |
| ERBB2 | Amplification may point to HER2-positive disease, but the label requires HER2 IHC 3+ | FDA label: Enhertu (fam-trastuzumab deruxtecan) 2026-05-15 | 2026-11-06 |
| FGFR1, FGFR2, FGFR3 | Activating mutation, fusion or rearrangement in tumor tissue | Recruiting trials open to any solid tumor after standard treatment: NCT06302621 (pemigatinib + afatinib), NCT06915753 (TYRA-430) | 2026-11-06 |

Labels retrieved from openFDA on 2026-11-06. One more tissue-agnostic indication can't be found by searching genes: tumor mutational burden-high (TMB-H, ≥10 mutations/Mb; Keytruda) needs a TMB measurement. If the record doesn't have one, it goes under Information gaps.

**Trials.** Recruiting tumor-agnostic trials are surveyed in [`tumor-agnostic-trials.md`](tumor-agnostic-trials.md) (2026-11-06: 34 genes with trials, e.g. KRAS G12 variants, MTAP loss, NRG1 fusions, PIK3CA H1047R, TP53 Y220C, CCNE1 amplification). All of them are in the OncoKB cancer gene list, so the driver sweep (§2.5) already searches them for every patient. The survey is used as a checklist instead: compare the patient's findings against it, and confirm each matching trial's required finding before listing it. A gene belongs in the table above only if it's outside the OncoKB list (so the sweep wouldn't search it) or has an FDA tissue-agnostic label.

### 2.4 Operations

Genomics server `https://fhir-gen-ops.herokuapp.com`; subject = the patient's MRN.

| Operation | Returns |
|---|---|
| `subject-operations/metadata-operations/$find-study-metadata` | The patient's tests: date, specimen, genome build, and the variant types each can detect. With `ranges`, also each test's studied and uncallable regions within those ranges |
| `utilities/get-feature-coordinates?gene=<symbol>` | The gene's range in GRCh37 and GRCh38, plus its transcripts |
| `utilities/normalize-variant?variant=<HGVS or SPDI>` | A variant's canonical form, for the specific-variant search |
| `subject-operations/genotype-operations/$find-subject-variants` | Simple variants (SNVs, indels) in a range |
| `subject-operations/genotype-operations/$find-subject-specific-variants` | Whether the patient has particular variants |
| `subject-operations/genotype-operations/$find-subject-structural-intersecting-variants` | Structural variants (deletions, duplications, inversions, copy-number changes) overlapping a range |
| `subject-operations/phenotype-operations/$find-subject-molecular-consequences` | Each variant's effect on the gene's protein |
| `subject-operations/phenotype-operations/$find-subject-tx-implications` | Therapeutic implications (CIViC, and the report's own) |
| `subject-operations/phenotype-operations/$find-subject-dx-implications` | Diagnostic implications, including ClinVar pathogenicity |

**Use the operations for the patient's variants, and the knowledge bases directly for what the variants mean.** The reference implementation holds only a small snapshot of CIViC and ClinVar, and molecular consequences only for some variants. So it's fine, and often better, to:

- query **CIViC** directly (civicdb.org API) for therapeutic implications
- query **ClinVar** directly (NCBI E-utilities) for pathogenicity and somatic oncogenicity
- compute **molecular consequences** directly (e.g. Ensembl VEP or SnpEff) for variants the server has none for

The operations remain the source for which variants the patient has, and for what was studied (`$find-study-metadata`). Cite the knowledge base queried and the date, as in §3.

**Failed liftover.** When a patient's tests are in different genome builds, the server lifts each query range over to the other build. Some ranges can't be lifted over, and the server then returns HTTP 422 "Failed LiftOver" for the whole request, so one bad range in a batch fails every gene in it. Per the Genomics Reporting IG's [liftover guidance](https://build.fhir.org/ig/HL7/genomics-reporting/operations.html#liftover), the client then adjusts the query region. Search the failing gene on its own; if it still fails, widen its range slightly, or split it into parts to isolate the segment that can't be lifted over and search the rest. Record any part left unsearched in the search log, and treat it as not tested.

**Search by gene, not by chromosome.** Request time grows with the number of variants returned (about 0.3 seconds per variant, tested 2026-11-06), and the server's host ends any request at 30 seconds with HTTP 503 while the query keeps running. For an exome-sized patient (HG00406), a 1-Mb range took 12 seconds and a 10-Mb range timed out, so whole-chromosome queries aren't practical. The server-side consequence filter (`$find-subject-molecular-consequences` with `featureConsequences`) cut the 1-Mb query to 5 seconds. Real exome or genome volumes would need a bulk operation on the server.

**Be gentle with the server.** It's a shared reference implementation. Send one request at a time, batch ranges (about 20 genes per request), pause between requests, and stop if it starts returning errors (e.g. HTTP 503).

### 2.5 Finding cancer drivers

A cancer driver is an alteration that helps drive the tumor's growth, whether or not any treatment targets it yet. Drivers are searched across all cancer genes, not just the minimum, and are summarized in the Cancer drivers section (§5). Most won't be druggable today. Recording them now means they can be rechecked as evidence grows (§8).

**Cancer genes:** `public/data/OncoKB_CancerGeneList.tsv`, the OncoKB Cancer Gene List (downloaded 2026-11-06 from https://www.oncokb.org/api/v1/utils/cancerGeneList.txt; free for academic research). Its `Gene Type` column gives each gene's role. Genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG are searched (858 genes in this download). Genes typed INSUFFICIENT_EVIDENCE or NEITHER, or with no type, are not.

**Driver rule.** A somatic variant in a cancer gene is a driver if it's classified oncogenic or likely oncogenic in ClinVar's somatic oncogenicity classifications, or is at a Cancer Hotspots hotspot, whatever the gene's role (e.g. TP53 R248Q, a missense hotspot in a tumor suppressor). Otherwise, it's a driver when the variant fits the gene's role:

| Gene role | Variant counts as a driver when it's | Determined with |
|---|---|---|
| Tumor suppressor (TSG) | Loss of function: truncating (stop-gained, frameshift), splice-site, or a structural variant that disrupts the gene (table below) | `$find-subject-molecular-consequences` (or VEP); structural variant search |
| Oncogene | Activating: amplification, a fusion, or a known activating deletion (table below) | Structural variant search |
| Both | Either of the above | |

**Structural variants.** An intersecting search returns every structural variant whose span covers the gene, so a single large event can appear under hundreds of genes. Reduce the results to distinct events first, then interpret each by where its breakpoints fall:

| Event | Tumor suppressor | Oncogene |
|---|---|---|
| Deletion inside the gene (both breakpoints in it) | Driver if it removes coding exons (check with VEP); not if intronic only | Driver only if it's a known activating deletion (e.g. EGFR exons 2–7, EGFRvIII; MET exon 14 splice region) |
| Deletion of the whole gene | Driver only if (a) copy number shows both copies lost, (b) the event is focal (under about 5 Mb), or (c) there's a second hit (e.g. a truncating variant) in the same gene. Otherwise report it as a possible single-copy loss, not a driver | Not a driver |
| Deletion with one breakpoint in the gene | Driver (the gene is broken) | Check the other breakpoint: if it's in another gene, a possible fusion (candidate driver) |
| Duplication inside the gene (both breakpoints in it) | Driver if it duplicates coding exons (it can disrupt the reading frame) | Not a driver unless a known activating duplication |
| Duplication covering the whole gene | Not a driver | Amplification only if copy number shows it (e.g. ≥ 5 copies, or as the trial or label defines it); otherwise not a driver |
| Inversion or translocation with a breakpoint in the gene | Driver (the gene is broken) | A possible fusion (candidate driver), especially if the other breakpoint is in a known fusion partner |
| Inversion spanning the gene with no breakpoint in it | Not a driver (the gene is moved, intact) | Not a driver |

The Genomics Explorer lists every structural variant that spans a gene, so genes merely spanned by a large event show rows there. When a reader would see such a gene in the Explorer, especially a gene they'd care about (e.g. one in the patient's germline panel or gene lists), say in Genomic findings that it's spanned, not broken, and why it isn't counted.

Breakpoints just outside a gene (in its promoter or regulatory region) aren't counted. Note them only if the literature describes that mechanism for the gene. If the test doesn't detect a variant type (e.g. no translocations in `dnaChangeType`), those events are not tested and go under Information gaps.

Other variants in cancer genes (for example a missense variant of unknown effect in a TSG) are not called drivers. Germline variants are not drivers in this sense; they stay in Genomic findings.

*Known limitation:* there's no comprehensive free source for activating variants in oncogenes. ClinVar's somatic classifications and Cancer Hotspots cover only part of them, so some activating variants may be missed. Where variant allele frequencies are missing or a placeholder (e.g. the same value for every variant), the synthesis can't say whether a driver is clonal (present in all tumor cells) or subclonal.

**Pathways:** where a driver's gene is in one of the TCGA pathway files (`public/data/pathways/`, Sanchez-Vega et al., *Cell* 2018), use that pathway to group it, as the Explorer's Pathways view does.

## 3. Evidence sources

**Look up, don't recall.** Background knowledge (a person's experience, or an AI model's training) may suggest what to look for: which genes matter, which drugs might apply, what resistance to expect. It is never the source of a fact in the synthesis. Every fact must come from a source looked up while writing this synthesis and cited with the date it was retrieved: the patient's data, label indications, evidence levels, trial status and eligibility, and literature. A training-based AI's knowledge also stops at a cutoff date, and anything newer reaches the synthesis only through a lookup. Anything that can't be verified goes under Information gaps rather than being stated from memory.

The exception is **general biology**: standard explanations of what a gene or pathway does (e.g. "PTEN loss activates PI3K signaling") may come from background knowledge, so the prose can explain without citing every sentence. Nothing specific to the patient (their variants, the variants' classification, evidence levels, or what applies to them) may come from memory.

| Source | What it's used for | How it's cited |
|---|---|---|
| FDA prescribing information (openFDA) | Whether a therapy is labeled for this patient's disease, setting and alteration | Brand name (generic name), label date, and the date the label was retrieved |
| CIViC | Therapeutic evidence for a variant, with its evidence level (A–E) and disease | Evidence item IDs (civicdb.org/evidence/&lt;id&gt;) |
| CIViC co-occurring-variant profiles (via `CIViC_Cooccurrence_KB.csv`) | Evidence for a combination of variants | Molecular profile ID, evidence item ID and level, and the date the knowledge base was retrieved |
| ClinVar | Whether a variant is pathogenic (germline) or oncogenic (somatic oncogenicity classifications) | "ClinVar records for &lt;variant&gt;" |
| OncoKB Cancer Gene List | Which genes are cancer genes, and whether each is an oncogene or tumor suppressor (§2.5) | OncoKB Cancer Gene List, and the date downloaded |
| Cancer Hotspots (cancerhotspots.org) | Whether an oncogene variant is a recurrent hotspot (§2.5) | Cancer Hotspots, and the date retrieved |
| TCGA pathways (Sanchez-Vega et al., *Cell* 2018) | Grouping drivers by pathway (§2.5) | The citation, as in `public/data/pathways/README.md` |
| ClinicalTrials.gov | Trials the patient might join | NCT IDs, the date searched, and recruitment status |
| Professional guidelines (e.g. ASCO) | Recommendations that aren't on a drug label, such as germline testing or hepatitis B management | Title, journal, year and DOI |
| Published literature (e.g. PubMed) and conference abstracts (e.g. ASCO, ESMO, AACR) | Mechanisms (such as resistance after a therapy), the genomic questions to query (§2.2), and new evidence not yet in labels or knowledge bases (the Emerging evidence section, §5). Not by itself a basis for listing a treatment option (see the rules below) | Title, journal or meeting, year and DOI or abstract number |
| FDA announcements (approvals, including accelerated approvals) | Approvals not yet reflected in the label retrieved from openFDA | Title and date of the announcement |

**Recent evidence.** Knowledge bases, labels and guidelines lag behind publication by months, and an AI's training stops at its cutoff date. So a definitive result published last month would be missed unless it's searched for. For each finding that matters (one that's actionable, or would be if evidence existed), and for the patient's disease setting, search for recent evidence:

- PubMed and the web, for the variant (or gene) in the patient's tumor type, limited to the search window
- FDA approval announcements
- abstracts from major oncology meetings (ASCO, ESMO, AACR)

The **search window** reaches back to the earliest of: the AI's training cutoff (if an AI is writing), the date `CIViC_Cooccurrence_KB.csv` and the CIViC data were retrieved, and the dates of the labels used. Record the window and the searches in the search log. What's found goes in Emerging evidence (§5). Searches can miss things, and some papers can be read only as abstracts, so the synthesis states the window searched rather than claiming nothing new exists.

**Rules for what is listed:**

- A therapy is listed as an option only when an FDA label covers the patient's situation, or a CIViC evidence item at level A or B does.
- Combination therapy for co-occurring alterations is an option only when an FDA label or CIViC co-occurring-variant evidence at level A or B supports it, and is then ranked with the other options (§5). Lower-level co-occurrence evidence can be mentioned, labeled with its level, in an option's considerations or under Not recommended (when it points to resistance), but not as a basis for treatment.
- New evidence not yet in a label or CIViC at level A or B (a recent trial result, an FDA announcement ahead of the label) goes in Emerging evidence, not in the treatment options, however strong it is.
- Evidence for a variant in a different disease is noted as such (e.g. "Evidence only in neuroblastoma"), not used for this patient.
- Only recruiting trials are offered as options. A trial the clinical record shows the patient is on is listed with its status on ClinicalTrials.gov, whatever it is.

## 4. Steps

1. **Gather the case.** Read the patient's bundle: any referral or visit reason, demographics, diagnosis and stage, biomarkers, treatment history and responses, progression, current status (ECOG, labs, comorbidities, allergies), and family history.
2. **Gather the genomic findings.** Take the reported variants and implications from the genetic reports. Check the patient's tests with `$find-study-metadata`, passing the ranges being searched. Search the genomics server for the minimum genes in §2.1, including the always-search list, and record each search in the search log (§2.2).
3. **Reconcile them.** Check that the two agree, and note where they don't: a reported variant not found in the genomic data, a finding in the genomic data not on the report, or a report implication the knowledge bases don't share. Keep germline and somatic findings apart, and note when origin isn't reported.
4. **Find cancer drivers.** Search the cancer genes (§2.5), apply the driver rule, and group the drivers by pathway. Drivers also go into the trial search (step 7).
5. **Look up evidence.** For each finding, record its CIViC evidence (level, disease, therapy) and ClinVar significance. Check the co-occurrence knowledge base for combinations of the patient's variants. Search for recent evidence over the search window (§3, Recent evidence).
6. **Look up treatment options.** For the patient's disease and treatment setting, check FDA labels for standard options, options matched to the genomic findings, and therapies that might seem to apply but don't (see §5, "Not indicated").
7. **Search for trials.** Search ClinicalTrials.gov for recruiting trials matching the genomic findings first, then the disease setting. Check the main eligibility criteria against the record. If a trial's criteria name a genomic finding that hasn't been searched, search it now (§2.2). If the trial is open to any tumor type, also add its gene to the always-search list (§2.3).
8. **Repeat steps 3–7** for any new finding, until the evidence raises no new genomic question (§2.2).
9. **Write** the synthesis using the template in §5 and the writing rules in §6.
10. **Check** it against §7 before saving.

## 5. Output template

Each synthesis has the same sections, in this order.

```markdown
# AI Case Synthesis: <Patient name> (MRN <MRN>)

> **Demo: synthesis pre-generated for this case** (<YYYY-MM-DD>). Research demonstration, synthetic patient, not for clinical use.

**Contents:** [Reason for visit](#reason-for-visit) · [Case summary](#case-summary) · [Genomic findings](#genomic-findings) · [Cancer drivers](#cancer-drivers) · [Information gaps](#information-gaps) · [Treatment options](#treatment-options) · [Emerging evidence](#emerging-evidence) · [References](#references)

## Reason for visit
## Case summary
## Genomic findings
## Cancer drivers
## Information gaps
## Treatment options
## Emerging evidence

---

## References
```

- **Contents:** one line of links to the sections, as above, under the banner. The AI Synthesis tab gives each section heading an id made from its text (lowercase, words joined by hyphens), so a renamed heading needs its link updated.
- **Reason for visit:** why the patient is being presented now, and the question or questions the tumor board is asked to address. One or two sentences on the trigger (e.g. a progression, a new biopsy or genetic report, a treatment decision point), then the questions as bullets. Keep it short: it points to the Case summary rather than repeating it. Infer the reason from the record as well as possible: the most recent events and the decision they lead to. If the record has a referral (e.g. a ServiceRequest for the tumor board, or an Encounter's reason), use its question.

  For example (L2345):

  > ## Reason for visit
  >
  > Progression on osimertinib (2026-09-28) after about 2.5 years of disease control; liver biopsy (2026-10-06) shows no transformation. New BRAF V600E, KRAS G12C and MET copy gain alongside persistent EGFR L858R/T790M.
  >
  > - What next-line therapy is supported?
  > - Should osimertinib be continued?
  > - Does the MET copy gain meet a trial's definition of amplification?
- **Case summary:** bullets, with bold lead-ins (Patient, Diagnosis, Treatment, Progression, Current status). Dates and values as recorded. Sub-bullets for lists of findings.
- **Genomic findings:** a table with columns Gene | Finding | Evidence, one row per finding. Bold the findings that carry evidence. A sentence after the table for anything the table can't say: how the findings fit together, how they changed over time, or that none has therapeutic evidence. When germline and tumor testing are both present, say which is which.
- **Cancer drivers:** a short section a clinician can take in at a glance, in this order:
  1. **Bottom line:** one bolded sentence on what drives the tumor and what that means for treatment.
  2. **Driver table:** Pathway | Gene | Role | Change | Treatment status. *Role* is "brake" (tumor suppressor, lost) or "accelerator" (oncogene, switched on). *Treatment status* is one of: targetable now (pointing to Treatment options), trial only, or watch (no drug yet; recheck as evidence grows). Major drivers first; minor ones (little known significance, no treatment link) last or left to the search log.
  3. **One paragraph of plain-language prose:** what's broken and what that does, how the drivers fit together and with the tumor type, and anything worth confirming (e.g. whether both copies of a gene are lost). Then say what can't be known (e.g. clonality).

  Writing for this section: use the brake and accelerator language consistently; name a pathway once, in plain words (e.g. "the PI3K growth pathway"); leave exon numbers, sizes and hotspot counts to Genomic findings and the search log; drivers already in Genomic findings are named, not described again. Mention the app's Pathways view for readers who want the pathway diagrams. If no drivers are found, say so in one sentence, with the genes searched.
- **Information gaps:** bullets for missing information that would change the options or the reading of the findings, saying why each matters.
- **Treatment options:** the options for the decision in Reason for visit, **ranked**, most preferred first. Each gives its rank, the label or evidence that supports it, then *Considerations* specific to this patient, and a short reason for its place in the ranking. Options for a later decision (e.g. after first-line treatment) follow under their own heading, ranked the same way, since they're a sequence and not alternatives to the current choice. The ranking weighs, in roughly this order:
  1. **Strength of evidence:** an FDA label for this disease and setting ahead of CIViC level A, ahead of level B; phase 3 evidence with a survival benefit ahead of single-arm or response-rate evidence. (Treatment guidelines will be added as a ranking source later.)
  2. **Fit to the patient's setting:** line of therapy, prior treatment (e.g. a drug already received or progressed on), and the label's or trial's eligibility.
  3. **Patient factors:** performance status, organ function, comorbidities, allergies and expected toxicity.

  Where options are too close to call, say so and rank them together rather than forcing an order. Clinical trials can be ranked among the options when one is a better fit than the standard choices; otherwise they're listed after the ranked options.

  **Combinations for co-occurring alterations are ranked with everything else**, not listed separately: one that addresses two of the patient's drivers can well rank first. Label such an option with the alterations it addresses (e.g. *"For co-occurring EGFR L858R + MET amplification"*) and cite the co-occurrence evidence. The same criteria apply, so strong co-occurrence evidence for a drug the patient has already progressed on ranks low or not at all. If the co-occurring-variant evidence was checked and none applies, say so in one line under Treatment options.

  Then, as they apply:
  - **Not indicated by the patient's findings:** therapies a reader might expect, and why the label or evidence doesn't cover this patient.
  - **Potential clinical trials:** a table (Trial | Design | Locations), headed with the search date and "eligibility must be confirmed with each site." If none matched, say what was searched.
- **Emerging evidence:** new evidence relevant to the patient's findings or disease setting that isn't yet in a label or CIViC at level A or B (§3, Recent evidence). Each item says what was found, its source and date, the kind of evidence (e.g. phase 3 trial, single-arm study, accelerated approval announcement), and how it bears on this patient, for the tumor board to weigh. It is not a recommendation. If nothing was found, say so and give the search window (e.g. "No new evidence found in searches covering 2026-06-01 to 2026-11-06").
- **References:** grouped by source (FDA labels, CIViC, co-occurrence, ClinVar, guidelines, literature, ClinicalTrials.gov), with the dates from §3. Also a **Genes searched** line, in plain words, summarizing what was searched and when (e.g. "small variants and copy-number changes in lung cancer genes and in 852 of 858 known cancer genes"), and anything that couldn't be searched. The details belong in the search log, which the synthesis doesn't name.

## 6. Writing rules

- Plain language for clinicians who aren't genomics specialists. Explain genomic terms they may not know (e.g. "acquired (somatic)").
- **Write only for the clinician.** The synthesis never mentions this README, the search log, file names, internal list names (e.g. "always-search list", "Actionable Gene List"), rule names ("driver rule") or section numbers. Describe what was done in clinical terms instead (e.g. "trials open to any solid tumor were matched against his alterations").
- Say when there is no evidence, explicitly, not by leaving a finding out.
- Keep germline and somatic findings apart. A somatic finding doesn't contradict a negative germline test, and the synthesis should say so where a reader might think it does.
- Tie every option to the patient: why they meet the label, or what would need to be confirmed (e.g. "a current AFP is needed").
- Recommending for or against an option is fine (e.g. "**Not recommended:** continuing osimertinib alone"), when the recommendation states the evidence behind it. Don't overstate the evidence: use "labeled for", "consistent with" and "would need" where they're accurate.
- Use the record's own values and dates. Don't fill in what the record doesn't say. Put it under Information gaps.

## 7. Checks before saving

- [ ] The genomic findings agree with the genetic reports and with what the Genomics Explorer shows for this patient, or the difference is explained.
- [ ] The minimum genes (§2.1), including the always-search list, were searched.
- [ ] The search log `<MRN>.searches.md` is saved, records every query with its reason, and the Genes searched line summarizes it.
- [ ] Every negative finding stated as absent was checked against the studied and uncallable regions and the variant types tested (§2.2). Anything not tested is under Information gaps.
- [ ] Every question in Reason for visit is answered in Treatment options or Emerging evidence, or left open under Information gaps.
- [ ] Every fact in the case summary is in the bundle.
- [ ] Every driver named in Cancer drivers meets the driver rule (§2.5), with its classification looked up; the section has a bottom line, a driver table and one paragraph of prose.
- [ ] Treatment options are ranked, each with a reason for its place.
- [ ] The co-occurring-variant evidence was checked for the patient's alterations, and Treatment options reflects it: qualifying combinations ranked and labeled, weaker evidence in considerations or Not recommended, or a line saying none applies.
- [ ] Every option, "not indicated" item and trial has a reference, and every reference is used.
- [ ] Every fact was looked up while writing this synthesis and is in the search log; none is stated from memory (§3).
- [ ] The recent-evidence search covered the search window (§3), and Emerging evidence states the window.
- [ ] Label, knowledge-base and trial dates are given, and trials were recruiting on the date searched.
- [ ] The banner line is present, with the date the synthesis was written.
- [ ] The file is named `<MRN>.md`, matching the MRN in the app's MRN box.

## 8. When to revise a synthesis

- The patient's bundle changes: new clinical data, a new genetic report, or new reported variants or implications.
- The genomic data or knowledge bases the app uses change in a way that affects the patient's findings.
- A drug label or a listed trial changes (e.g. a trial stops recruiting).
- This README's process changes.

When revising, update the banner date and the retrieval dates of anything looked up again.
