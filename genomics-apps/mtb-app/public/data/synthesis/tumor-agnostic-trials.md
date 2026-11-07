# Tumor-agnostic trial survey

Recruiting ClinicalTrials.gov studies (API v2), with at least one US site, whose condition or title is generic (e.g. "solid tumor", "advanced cancer") and whose inclusion criteria name the gene. Surveyed 2026-11-06 for the always-search list (`README.md` §2.3).

This is a screen, not a verified list: matches are by text, and some are disease-specific cohorts within a broader study (e.g. EGFR in NSCLC). When a patient has a finding in one of these genes, open the trials and confirm the required finding (variant class, origin, prior therapy) before listing one.

All genes below are in the OncoKB cancer gene list, so the cancer-driver sweep (§2.5) already searches them for every patient.

| Gene | Trials matched | Examples (inclusion text) |
|---|---|---|
| KRAS | 44 | NCT06917079: "KRAS G12A, G12C, G12D, G12S, or G12V mutation"; NCT07397338: "KRAS, HRAS, or NRAS" |
| EGFR | 18 | NCT03175224: "EGFR NSCLC: NSCLC regardless of histology, harboring MET amplification and wild-type"; NCT07444814: "EGFR-Wt NSCLC" |
| MTAP | 15 | NCT06630416: "MTAP loss of function mutation including but not limited to SNVs that cause truncati"; NCT07382544: "MTAP as detected by Clinical Laboratory Improvement Amendmentscertified next-generat" |
| NRAS | 8 | NCT07397338: "NRAS"; NCT06026410: "NRAS, and/or HRAS-mutant and/or amplified NSCLC or CRC; KRAS-mutant and/or amplified" |
| BRAF | 7 | NCT05503797: "BRAF gene fusion in tumor and/or blood detected by an analytically validated test by"; NCT06326411: "BRAF alterations" |
| PIK3CA | 6 | NCT07819136: "PIK3CA H1047R mutation in tumor tissue and/or blood (circulating tumor deoxyribonuclei"; NCT07524322: "PIK3CA mutation in tumor tissue and/or blood" |
| TP53 | 6 | NCT04585750: "TP53 Y220C mutation"; NCT07060989: "TP53 Y220C mutation" |
| ROS1 | 5 | NCT04094610: "ROS1 point mutation, fusion, or amplification or NTRK1-3 fusion as identified by loc"; NCT03093116: "ROS1, NTRK1, NTRK2, or NTRK3 gene rearrangement by protocol specified tests" |
| ALK | 4 | NCT05384626: "ALK rearrangement or activating ALK mutation"; NCT03093116: "ALK, ROS1, NTRK1, NTRK2, or NTRK3 gene rearrangement by protocol specified tests" |
| CCNE1 | 4 | NCT06586957: "CCNE1 amplification)"; NCT05867251: "CCNE1-amplified tumors that are either epithelial ovarian cancer, primary peritoneal" |
| HRAS | 4 | NCT07397338: "HRAS, or NRAS"; NCT06026410: "HRAS-mutant and/or amplified tumors (any solid tumor type); HRAS overexpression (onl" |
| ATM | 3 | NCT04657068: "ATM gene likely to predict for loss of ATM protein"; NCT05694715: "ATM, and/or PALB2 (based upon archival tumor tissue or germ line testing from any C" |
| BRCA1 | 3 | NCT07356453: "BRCA1/2 mutations (germline or somatic, pathogenic or likely pathogenic) and other mu"; NCT04890613: "BRCA1/2 and/or PALB2, and/or any genetic alterations listed below as indicated in a C" |
| BRCA2 | 3 | NCT04890613: "BRCA2, BRIP1, CCNE1, CHEK1, CHEK2, CDK12, CREBBP, FANCA, FANCI, FANCL, FANC2, FANCB,"; NCT05694715: "BRCA2, ATM, and/or PALB2 (based upon archival tumor tissue or germ line testing from" |
| ERBB2 | 3 | NCT05673928: "ERBB2 mutation(s) (verified by MDACC Precision Oncology Decision Support)"; NCT07361562: "ERBB2-activating alteration or NRG1 gene fusion in blood and/or tumor or HER2 overexp" |
| MET | 3 | NCT06669117: "MET alterations based on local documentation of blood or archived tissue results:"; NCT03175224: "MET inhibitor naive in first line: Histologically or cytologically confirmed NSCLC" |
| NF1 | 3 | NCT03454035: "NF1 mutations, subjects with loss-of-function NF1mutations and without any BRAFV600"; NCT02390752: "NF1 and with malignant peripheral nerve sheath tumor (MPNST)" |
| FGFR2 | 2 | NCT07359820: "FGFR2 gene fusion or rearrangement per local testing of blood and/or tumor"; NCT06302621: "FGFR2 extracellular domain in-frame deletions on tumor profiling in tumor tissue as d" |
| NTRK1 | 2 | NCT03093116: "NTRK1, NTRK2, or NTRK3 gene rearrangement by protocol specified tests"; NCT03556228: "NTRK1 gene fusion positive ("NTRK1+") solid tumors or lymphomas, that is relapsed, re" |
| PALB2 | 2 | NCT04890613: "PALB2, and/or any genetic alterations listed below as indicated in a Clinical Laborat"; NCT05694715: "PALB2 (based upon archival tumor tissue or germ line testing from any Clinical Labora" |
| PTEN | 2 | NCT04997993: "PTEN expression as determined by immunohistochemistry"; NCT07109726: "PTEN pathway alteration detected by a sponsor approved test" |
| SMARCA4 | 2 | NCT06561685: "SMARCA4 (BRG1) alteration:"; NCT07284186: "SMARCA4 loss-of-function mutation that have progressed on, are intolerant of, or not ot" |
| ARID1A | 1 | NCT04104776: "ARID1A mutant urothelial carcinoma or other ARID1A mutant solid tumors (with cohort sp" |
| BAP1 | 1 | NCT04104776: "BAP1 loss" |
| FGFR3 | 1 | NCT06874335: "FGFR3 alteration (mutation or fusion)" |
| FGFR4 | 1 | NCT06915753: "FGFR4 gene, or focal amplifications of FGF19" |
| IDH1 | 1 | NCT05417594: "IDH1/2-mutant glioma" |
| KMT2D | 1 | NCT06630416: "KMT2D-gene) or UTX-protein (KDM6A-gene) or MTAP loss of function mutation including b" |
| MDM2 | 1 | NCT06578624: "MDM2 amplification); OR Cohort B: Tumor p53 wild-type in other solid tumor" |
| NF2 | 1 | NCT06251310: "NF2 mutations" |
| NRG1 | 1 | NCT07361562: "NRG1 gene fusion in blood and/or tumor or HER2 overexpression in tumor" |
| RET | 1 | NCT04683250: "RET-gene abnormalities determined on tissue biopsy or liquid biopsy" |
| STK11 | 1 | NCT05098132: "STK11 mutated on tumor tissue or ctDNA" |
| XPO1 | 1 | NCT05985161: "XPO1 Dependency: Defined as either Darwin OncoTarget demonstrating XPO1 as aberrantl" |

Searched with no matching trial: CDKN2A, ERBB3, AKT1, IDH2, ATR, SMARCB1, KEAP1, FGFR1, MYC, CDK12, POLE, TSC1, TSC2, PBRM1.

Survey method: for each gene, `query.term=<gene>`, `query.cond=solid tumor`, `filter.overallStatus=RECRUITING`; kept studies with a US site, a generic condition or title, and the gene in the inclusion criteria.
