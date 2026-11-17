# Search log: TCGA-DD-A1EH

Companion to the case synthesis. Every lookup made while writing it on 2026-11-16, with why and the result. Genomics server: https://fhir-gen-ops.herokuapp.com, subject TCGA-DD-A1EH. Requests were sent one at a time with 2 s pauses; another patient's synthesis was being written at the same time.

## Clinical record

Patient transaction bundle (clinical-data/TCGA-DD-A1EH/TCGA-DD-A1EH-transaction-bundle.json), all 57 entries read: Patient, Conditions (HCC, liver recurrence, lung metastases, chronic hepatitis B), TNM, Procedures (liver biopsy 2025-12-20, right hepatectomy 2026-01-26), DiagnosticReports (CT 2026-01-06, pathology 2026-01-26, CT 2026-03-11, CT 2026-05-02, CT 2026-07-13, MRI 2026-07-14, genomics report 2026-07-28), Observations (AFP, Child-Pugh, ECOG, MSI, TMB, disease status, variants, therapeutic implication), MedicationRequests (tenofovir, sorafenib, tazemetostat), ResearchSubject/ResearchStudy (NCT05023655), AllergyIntolerance (latex), ImagingStudies. No referral or encounter reason in the record; reason for visit inferred from the trial and drug status.

## Genome queries

### Round 1: minimum genes and cancer-gene sweep

| # | Operation | Query | Why | Result |
|---|---|---|---|---|
| 1 | `$find-study-metadata` (no ranges) | subject only | Tests on file | One test TCGA-DD-A1EH-T1, date 2020-01-20, specimen TCGA-DD-A1EH-Sp1, GRCh37; dnaChangeType SNV, MNV, delins (no CNV/SV types). Note: insertions are nevertheless returned (ARID1A, AXIN1), so the listed types understate what the test holds |
| – | `utilities/get-feature-coordinates?gene=` × 869 | Every gene below, one request each | GRCh37 ranges for the searches | 863 resolved; no coordinates for RYBP, NHERF1, TRA, TRB, TRD, TRG (resolved in round 2) |
| 2 | `$find-subject-variants` | ARID1A, ARID2, AXIN1, CD44, CDK9, CEACAM1, DNAJB1, EGFR, FGF19, FGF3, FGF4, FGFR1, FGFR2, FGFR4, GALNT14, JAK1, KIF23, KRAS, NEDD9, PBLD (GRCh37 gene ranges) | Minimum: reported genes (ARID1A, JAK1), HCC Extended Gene List (25 genes), tumor-agnostic always-search genes (NTRK1-3, BRAF, RET, MLH1, MSH2, MSH6, PMS2, ERBB2, FGFR1-3; FGFR rows since removed from the list but searched as cancer genes), plus the first OncoKB genes | 4 variant(s) |
| 3 | `$find-subject-structural-intersecting-variants` | ARID1A, ARID2, AXIN1, CD44, CDK9, CEACAM1, DNAJB1, EGFR, FGF19, FGF3, FGF4, FGFR1, FGFR2, FGFR4, GALNT14, JAK1, KIF23, KRAS, NEDD9, PBLD (GRCh37 gene ranges) | Minimum: reported genes (ARID1A, JAK1), HCC Extended Gene List (25 genes), tumor-agnostic always-search genes (NTRK1-3, BRAF, RET, MLH1, MSH2, MSH6, PMS2, ERBB2, FGFR1-3; FGFR rows since removed from the list but searched as cancer genes), plus the first OncoKB genes | 0 variant(s) |
| 4 | `$find-subject-variants` | PRKACA, PROM1, PTPN6, RAF1, VEGFC, NTRK1, NTRK2, NTRK3, BRAF, RET, MLH1, MSH2, MSH6, PMS2, ERBB2, FGFR3, ABL1, AKT1, ALK, AMER1 (GRCh37 gene ranges) | Minimum: reported genes (ARID1A, JAK1), HCC Extended Gene List (25 genes), tumor-agnostic always-search genes (NTRK1-3, BRAF, RET, MLH1, MSH2, MSH6, PMS2, ERBB2, FGFR1-3; FGFR rows since removed from the list but searched as cancer genes), plus the first OncoKB genes | 0 variant(s) |
| 5 | `$find-subject-structural-intersecting-variants` | PRKACA, PROM1, PTPN6, RAF1, VEGFC, NTRK1, NTRK2, NTRK3, BRAF, RET, MLH1, MSH2, MSH6, PMS2, ERBB2, FGFR3, ABL1, AKT1, ALK, AMER1 (GRCh37 gene ranges) | Minimum: reported genes (ARID1A, JAK1), HCC Extended Gene List (25 genes), tumor-agnostic always-search genes (NTRK1-3, BRAF, RET, MLH1, MSH2, MSH6, PMS2, ERBB2, FGFR1-3; FGFR rows since removed from the list but searched as cancer genes), plus the first OncoKB genes | 0 variant(s) |
| 6 | `$find-subject-variants` | APC, AR, ASXL1, ATM, ATRX, BAP1, BCL2, BCOR, BRCA1, BRCA2, CARD11, CBL, CDC73, CDH1, CDKN2A, CEBPA, CIC, CREBBP, CTNNB1, DAXX (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 7 | `$find-subject-structural-intersecting-variants` | APC, AR, ASXL1, ATM, ATRX, BAP1, BCL2, BCOR, BRCA1, BRCA2, CARD11, CBL, CDC73, CDH1, CDKN2A, CEBPA, CIC, CREBBP, CTNNB1, DAXX (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 8 | `$find-subject-variants` | DNMT3A, EP300, EZH2, FBXW7, FLT3, FOXL2, GATA3, GNA11, GNAQ, GNAS, HNF1A, HRAS, IDH1, IDH2, JAK2, JAK3, KDM5C, KDM6A, KIT, KMT2D (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 9 | `$find-subject-structural-intersecting-variants` | DNMT3A, EP300, EZH2, FBXW7, FLT3, FOXL2, GATA3, GNA11, GNAQ, GNAS, HNF1A, HRAS, IDH1, IDH2, JAK2, JAK3, KDM5C, KDM6A, KIT, KMT2D (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 10 | `$find-subject-variants` | MAP2K1, MAP3K1, MED12, MEN1, MET, MPL, MYD88, NF1, NF2, NFE2L2, NOTCH1, NOTCH2, NPM1, NRAS, PAX5, PBRM1, PDGFRA, PIK3CA, PIK3R1, PPP2R1A (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 11 | `$find-subject-structural-intersecting-variants` | MAP2K1, MAP3K1, MED12, MEN1, MET, MPL, MYD88, NF1, NF2, NFE2L2, NOTCH1, NOTCH2, NPM1, NRAS, PAX5, PBRM1, PDGFRA, PIK3CA, PIK3R1, PPP2R1A (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 12 | `$find-subject-variants` | PRDM1, PTCH1, PTEN, PTPN11, RB1, RNF43, SETD2, SF3B1, SMAD2, SMAD4, SMARCA4, SMARCB1, SMO, SOCS1, SPOP, STAG2, STK11, TET2, TNFAIP3, TP53 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 13 | `$find-subject-structural-intersecting-variants` | PRDM1, PTCH1, PTEN, PTPN11, RB1, RNF43, SETD2, SF3B1, SMAD2, SMAD4, SMARCA4, SMARCB1, SMO, SOCS1, SPOP, STAG2, STK11, TET2, TNFAIP3, TP53 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 14 | `$find-subject-variants` | TSC1, U2AF1, VHL, WT1, AKT2, ATR, B2M, BARD1, BCL6, BRD4, BRIP1, BTK, CALR, CASP8, CBFB, CCND1, CCND2, CCND3, CCNE1, CD274 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 15 | `$find-subject-structural-intersecting-variants` | TSC1, U2AF1, VHL, WT1, AKT2, ATR, B2M, BARD1, BCL6, BRD4, BRIP1, BTK, CALR, CASP8, CBFB, CCND1, CCND2, CCND3, CCNE1, CD274 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 16 | `$find-subject-variants` | CD79A, CD79B, CDK12, CDK4, CDK6, CDKN1B, CDKN2C, CHEK2, CRLF2, CSF3R, CTCF, CXCR4, DDR2, ERBB3, ERBB4, ERG, ESR1, ETV6, FANCA, FANCC (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 17 | `$find-subject-structural-intersecting-variants` | CD79A, CD79B, CDK12, CDK4, CDK6, CDKN1B, CDKN2C, CHEK2, CRLF2, CSF3R, CTCF, CXCR4, DDR2, ERBB3, ERBB4, ERG, ESR1, ETV6, FANCA, FANCC (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 18 | `$find-subject-variants` | FLCN, FUBP1, GATA2, H3-3A, H3C2, HGF, IRF4, JUN, KDM5A, KDR, KEAP1, KMT2A, KMT2C, MAP2K2, MAP2K4, MAPK1, MDM2, MDM4, MITF, MTOR (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 19 | `$find-subject-structural-intersecting-variants` | FLCN, FUBP1, GATA2, H3-3A, H3C2, HGF, IRF4, JUN, KDM5A, KDR, KEAP1, KMT2A, KMT2C, MAP2K2, MAP2K4, MAPK1, MDM2, MDM4, MITF, MTOR (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 20 | `$find-subject-variants` | MUTYH, MYC, MYCL, MYCN, NKX2-1, NSD2, PALB2, PDCD1LG2, PDGFRB, PHF6, PRKAR1A, RAD21, ROS1, RUNX1, SDHA, SDHB, SDHC, SDHD, SOX2, SPEN (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 21 | `$find-subject-structural-intersecting-variants` | MUTYH, MYC, MYCL, MYCN, NKX2-1, NSD2, PALB2, PDCD1LG2, PDGFRB, PHF6, PRKAR1A, RAD21, ROS1, RUNX1, SDHA, SDHB, SDHC, SDHD, SOX2, SPEN (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 22 | `$find-subject-variants` | SRC, STAT3, SUFU, SYK, TENT5C, TGFBR2, TNFRSF14, TSC2, XPO1, AKT3, ARAF, ARID1B, AURKA, AURKB, AXL, BCL10, BCORL1, BCR, BIRC3, BLM (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 23 | `$find-subject-structural-intersecting-variants` | SRC, STAT3, SUFU, SYK, TENT5C, TGFBR2, TNFRSF14, TSC2, XPO1, AKT3, ARAF, ARID1B, AURKA, AURKB, AXL, BCL10, BCORL1, BCR, BIRC3, BLM (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 24 | `$find-subject-variants` | BTG1, CDK8, CDKN2B, CHEK1, CRKL, CYLD, DOT1L, EED, EIF4A2, EPHA3, EPHB1, ERCC4, ETV1, FAS, FH, FLT1, FLT4, FOXO1, FOXP1, GRIN2A (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 25 | `$find-subject-structural-intersecting-variants` | BTG1, CDK8, CDKN2B, CHEK1, CRKL, CYLD, DOT1L, EED, EIF4A2, EPHA3, EPHB1, ERCC4, ETV1, FAS, FH, FLT1, FLT4, FOXO1, FOXP1, GRIN2A (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 26 | `$find-subject-variants` | GSK3B, ID3, IGF1R, IKBKE, IL7R, INPP4B, IRS2, KLF4, LMO1, MAP3K13, MCL1, MEF2B, MRE11, MSH3, MSI2, NBN, NCOR1, NFKBIA, NSD1, NT5C2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 27 | `$find-subject-structural-intersecting-variants` | GSK3B, ID3, IGF1R, IKBKE, IL7R, INPP4B, IRS2, KLF4, LMO1, MAP3K13, MCL1, MEF2B, MRE11, MSH3, MSI2, NBN, NCOR1, NFKBIA, NSD1, NT5C2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 28 | `$find-subject-variants` | P2RY8, PDCD1, PIK3CB, POLD1, POLE, POT1, PPARG, RAC1, RAD50, RAD51, RAD51B, RBM10, REL, RHOA, RICTOR, RPTOR, SETBP1, SOX9, STAT5B, SUZ12 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 29 | `$find-subject-structural-intersecting-variants` | P2RY8, PDCD1, PIK3CB, POLD1, POLE, POT1, PPARG, RAC1, RAD50, RAD51, RAD51B, RBM10, REL, RHOA, RICTOR, RPTOR, SETBP1, SOX9, STAT5B, SUZ12 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 30 | `$find-subject-variants` | TBX3, TCF3, TERT, TET1, TOP1, TP63, TSHR, ACVR1, ACVR1B, ALOX12B, AXIN2, BCL11B, BMPR1A, CDKN1A, CIITA, CUL3, CUX1, DDX3X, DICER1, DIS3 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 31 | `$find-subject-structural-intersecting-variants` | TBX3, TCF3, TERT, TET1, TOP1, TP63, TSHR, ACVR1, ACVR1B, ALOX12B, AXIN2, BCL11B, BMPR1A, CDKN1A, CIITA, CUL3, CUX1, DDX3X, DICER1, DIS3 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 32 | `$find-subject-variants` | DNMT1, EPHA7, ERCC2, ERCC3, ERCC5, ERRFI1, ETV4, ETV5, EWSR1, FANCD2, FANCG, FAT1, FBXO11, FOXA1, GLI1, GNA13, HDAC1, HLA-A, IKZF3, INHBA (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 1 variant(s) |
| 33 | `$find-subject-structural-intersecting-variants` | DNMT1, EPHA7, ERCC2, ERCC3, ERCC5, ERRFI1, ETV4, ETV5, EWSR1, FANCD2, FANCG, FAT1, FBXO11, FOXA1, GLI1, GNA13, HDAC1, HLA-A, IKZF3, INHBA (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 34 | `$find-subject-variants` | LATS1, LATS2, LYN, MAP3K14, MAX, MST1R, MYB, MYOD1, NOTCH3, PARP1, PHOX2B, PIK3CG, PIK3R2, PLCG2, PPM1D, PPP6C, PRKCI, PRKN, PTPRT, RAD51C (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 35 | `$find-subject-structural-intersecting-variants` | LATS1, LATS2, LYN, MAP3K14, MAX, MST1R, MYB, MYOD1, NOTCH3, PARP1, PHOX2B, PIK3CG, PIK3R2, PLCG2, PPM1D, PPP6C, PRKCI, PRKN, PTPRT, RAD51C (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 36 | `$find-subject-variants` | RAD51D, RECQL4, RRAS2, RUNX1T1, SDHAF2, SGK1, SH2B3, SMAD3, STAT5A, STAT6, TBL1XR1, TCF7L2, TMEM127, VEGFA, WWTR1, XRCC2, ZFHX3, ABI1, ABL2, ABRAXAS1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 37 | `$find-subject-structural-intersecting-variants` | RAD51D, RECQL4, RRAS2, RUNX1T1, SDHAF2, SGK1, SH2B3, SMAD3, STAT5A, STAT6, TBL1XR1, TCF7L2, TMEM127, VEGFA, WWTR1, XRCC2, ZFHX3, ABI1, ABL2, ABRAXAS1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 38 | `$find-subject-variants` | AFDN, AFF4, ANKRD11, ARHGAP26, ARHGAP35, ARHGEF12, ARID5B, ARNT, ASXL2, ATF1, ATIC, BBC3, BCL11A, BCL2L11, BCL2L2, BCL3, BCL7A, BCL9, BTG2, CAMTA1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 39 | `$find-subject-structural-intersecting-variants` | AFDN, AFF4, ANKRD11, ARHGAP26, ARHGAP35, ARHGEF12, ARID5B, ARNT, ASXL2, ATF1, ATIC, BBC3, BCL11A, BCL2L11, BCL2L2, BCL3, BCL7A, BCL9, BTG2, CAMTA1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 40 | `$find-subject-variants` | CBFA2T3, CCNQ, CD276, CD58, CD70, CD74, CDC42, CDX2, CEP43, COL1A1, COP1, CTLA4, CTNNA1, CYSLTR2, DCUN1D1, DEK, DNM2, DNMT3B, DTX1, DUSP22 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 41 | `$find-subject-structural-intersecting-variants` | CBFA2T3, CCNQ, CD276, CD58, CD70, CD74, CDC42, CDX2, CEP43, COL1A1, COP1, CTLA4, CTNNA1, CYSLTR2, DCUN1D1, DEK, DNM2, DNMT3B, DTX1, DUSP22 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 42 | `$find-subject-variants` | DUSP4, E2F3, EGFL7, EIF4E, ELF3, ELF4, ELL, EML4, EMSY, EPCAM, ERF, EZH1, EZR, FANCE, FANCF, FANCL, FCGR2B, FGF10, FGF14, FGF23 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 43 | `$find-subject-structural-intersecting-variants` | DUSP4, E2F3, EGFL7, EIF4E, ELF3, ELF4, ELL, EML4, EMSY, EPCAM, ERF, EZH1, EZR, FANCE, FANCF, FANCL, FCGR2B, FGF10, FGF14, FGF23 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 44 | `$find-subject-variants` | FGF6, FHIT, FLI1, FUS, FYN, GNA12, GNB1, GPS2, GREM1, H1-3, H1-4, H2AC17, H3C3, HDAC4, HDAC7, HIF1A, HIP1, HLA-B, HMGA1, HMGA2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 45 | `$find-subject-structural-intersecting-variants` | FGF6, FHIT, FLI1, FUS, FYN, GNA12, GNB1, GPS2, GREM1, H1-3, H1-4, H2AC17, H3C3, HDAC4, HDAC7, HIF1A, HIP1, HLA-B, HMGA1, HMGA2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 46 | `$find-subject-variants` | HOXA11, HOXA13, HOXB13, HSP90AA1, ICOSLG, IFNGR1, IGF1, IGF2, INHA, INPPL1, INSR, IRF1, IRF8, IRS1, JARID2, KLHL6, KMT2B, LASP1, LCK, LEF1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 47 | `$find-subject-structural-intersecting-variants` | HOXA11, HOXA13, HOXB13, HSP90AA1, ICOSLG, IFNGR1, IGF1, IGF2, INHA, INPPL1, INSR, IRF1, IRF8, IRS1, JARID2, KLHL6, KMT2B, LASP1, LCK, LEF1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 48 | `$find-subject-variants` | LMO2, LRP1B, LZTR1, MAFB, MAPK3, MECOM, MGA, MLF1, MLLT10, MLLT3, MN1, MSI1, MST1, MTAP, MYH11, NCOA3, NCSTN, NKX3-1, NOTCH4, NR4A3 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 49 | `$find-subject-structural-intersecting-variants` | LMO2, LRP1B, LZTR1, MAFB, MAPK3, MECOM, MGA, MLF1, MLLT10, MLLT3, MN1, MSI1, MST1, MTAP, MYH11, NCOA3, NCSTN, NKX3-1, NOTCH4, NR4A3 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 50 | `$find-subject-variants` | NTHL1, NUF2, NUP214, NUP98, PAFAH1B2, PAK1, PAK5, PAX3, PAX7, PBX1, PDE4DIP, PDGFB, PDK1, PGR, PICALM, PIK3CD, PIK3R3, PLCG1, PMAIP1, PML (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 51 | `$find-subject-structural-intersecting-variants` | NTHL1, NUF2, NUP214, NUP98, PAFAH1B2, PAK1, PAK5, PAX3, PAX7, PBX1, PDE4DIP, PDGFB, PDK1, PGR, PICALM, PIK3CD, PIK3R3, PLCG1, PMAIP1, PML (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 52 | `$find-subject-variants` | PMS1, PTPN2, PTPRD, PTPRS, QKI, RAB35, RAC2, RASA1, RBM15, RECQL, RHEB, RIT1, RPS6KA4, RPS6KB2, RRAGC, RRAS, RSPO2, RTEL1, SESN1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 53 | `$find-subject-structural-intersecting-variants` | PMS1, PTPN2, PTPRD, PTPRS, QKI, RAB35, RAC2, RASA1, RBM15, RECQL, RHEB, RIT1, RPS6KA4, RPS6KB2, RRAGC, RRAS, RSPO2, RTEL1, SESN1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 54 | `$find-subject-variants` | SESN2, SESN3, SET, SETDB1, SH2D1A, SHOC2, SHQ1, SLX4, SMC1A, SMC3, SMYD3, SOS1, SOX17, SPRED1, SS18, STK19, TAL1, TCL1A, TFE3, TGFBR1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 55 | `$find-subject-structural-intersecting-variants` | SESN2, SESN3, SET, SETDB1, SH2D1A, SHOC2, SHQ1, SLX4, SMC1A, SMC3, SMYD3, SOS1, SOX17, SPRED1, SS18, STK19, TAL1, TCL1A, TFE3, TGFBR1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 56 | `$find-subject-variants` | TLX1, TLX3, TNFRSF17, TP53BP1, TPM3, TRAF3, TRAF5, TYK2, UBR5, USP6, USP8, XBP1, XIAP, YAP1, YES1, ZNF217, ZNF384, ACKR3, ACSL3, ACTG1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 57 | `$find-subject-structural-intersecting-variants` | TLX1, TLX3, TNFRSF17, TP53BP1, TPM3, TRAF3, TRAF5, TYK2, UBR5, USP6, USP8, XBP1, XIAP, YAP1, YES1, ZNF217, ZNF384, ACKR3, ACSL3, ACTG1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 58 | `$find-subject-variants` | ACVR2A, ADGRA2, AGO1, APLNR, APOBEC3B, ARHGEF28, ARID3A, ARID3B, ARID4A, ARID4B, ATG5, ATP1A1, ATP6V1B2, ATXN2, ATXN7, BACH2, BCL9L, BRD3, BRSK1, BUB1B (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 59 | `$find-subject-structural-intersecting-variants` | ACVR2A, ADGRA2, AGO1, APLNR, APOBEC3B, ARHGEF28, ARID3A, ARID3B, ARID4A, ARID4B, ATG5, ATP1A1, ATP6V1B2, ATXN2, ATXN7, BACH2, BCL9L, BRD3, BRSK1, BUB1B (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 60 | `$find-subject-variants` | CACNA1D, CAD, CANT1, CBLB, CCN6, CD28, CDH11, CHD2, CHD4, CILK1, CLIP1, CMTR2, COL2A1, CRBN, CREB1, CRTC1, CTR9, CUL4A, CYP19A1, DDB2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 61 | `$find-subject-structural-intersecting-variants` | CACNA1D, CAD, CANT1, CBLB, CCN6, CD28, CDH11, CHD2, CHD4, CILK1, CLIP1, CMTR2, COL2A1, CRBN, CREB1, CRTC1, CTR9, CUL4A, CYP19A1, DDB2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 62 | `$find-subject-variants` | DDX5, DUSP9, ECT2L, EGR1, ELK4, ELN, EP400, EPHB4, EPOR, ERC1, ESCO2, ETAA1, ETS1, EXT1, EZHIP, FES, FOXF1, FRS2, FURIN, GAB1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 63 | `$find-subject-structural-intersecting-variants` | DDX5, DUSP9, ECT2L, EGR1, ELK4, ELN, EP400, EPHB4, EPOR, ERC1, ESCO2, ETAA1, ETS1, EXT1, EZHIP, FES, FOXF1, FRS2, FURIN, GAB1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 64 | `$find-subject-variants` | GAB2, GABRA6, GATA4, GATA6, GPC3, GTF2I, H1-5, HIRA, HLA-C, HSD3B1, IL3, IL6ST, IRF2, JAZF1, KEL, KLF5, KSR2, LMNA, LTB, LTK (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 1 variant(s) |
| 65 | `$find-subject-structural-intersecting-variants` | GAB2, GABRA6, GATA4, GATA6, GPC3, GTF2I, H1-5, HIRA, HLA-C, HSD3B1, IL3, IL6ST, IRF2, JAZF1, KEL, KLF5, KSR2, LMNA, LTB, LTK (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 66 | `$find-subject-variants` | MAD2L2, MAP3K7, MEF2C, MERTK, MGAM, MKI67, MNX1, MOB3B, MYO5A, NAB2, NADK, NCOA4, NFATC2, NFE2, NFKBIE, NRG1, PAK3, PAX8, PDS5B, PGBD5 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 67 | `$find-subject-structural-intersecting-variants` | MAD2L2, MAP3K7, MEF2C, MERTK, MGAM, MKI67, MNX1, MOB3B, MYO5A, NAB2, NADK, NCOA4, NFATC2, NFE2, NFKBIE, NRG1, PAK3, PAX8, PDS5B, PGBD5 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 68 | `$find-subject-variants` | PIGA, PIK3C2B, PPP2R2A, PRCC, PRKDC, PTK7, PTPN1, PTPN13, RANBP2, REST, ROBO1, RPL5, RSPO3, SAMHD1, SCG5, SETD1A, SETDB2, SLFN11, SMARCA1, SMARCA2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 69 | `$find-subject-structural-intersecting-variants` | PIGA, PIK3C2B, PPP2R2A, PRCC, PRKDC, PTK7, PTPN1, PTPN13, RANBP2, REST, ROBO1, RPL5, RSPO3, SAMHD1, SCG5, SETD1A, SETDB2, SLFN11, SMARCA1, SMARCA2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 70 | `$find-subject-variants` | SMG1, SOCS3, SOX10, SP140, SPRTN, STAG1, STAT4, TAF1, TET3, TRIM27, TRIP13, VAV1, VAV2, WIF1, WRN, XPA (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 71 | `$find-subject-structural-intersecting-variants` | SMG1, SOCS3, SOX10, SP140, SPRTN, STAG1, STAT4, TAF1, TET3, TRIM27, TRIP13, VAV1, VAV2, WIF1, WRN, XPA (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 72 | `$find-subject-variants` | XPC, YWHAE, ZBTB7A, ZMYM3, ZNRF3, A1CF, AAMP, ABCB1, ABCC3, ACACA, ACP3, ACSL4, ADARB2, ADHFE1, AGGF1, AGK, AIM2, AIP, AJUBA, AKT1S1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 73 | `$find-subject-structural-intersecting-variants` | XPC, YWHAE, ZBTB7A, ZMYM3, ZNRF3, A1CF, AAMP, ABCB1, ABCC3, ACACA, ACP3, ACSL4, ADARB2, ADHFE1, AGGF1, AGK, AIM2, AIP, AJUBA, AKT1S1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 74 | `$find-subject-variants` | ALDH1A3, ALDH1L2, ALDH2, ALOX15B, ALOX5, ANK1, ANKRD26, APCDD1, APEX1, ART1, ASCL1, ASXL3, ATMIN, ATRIP, AURKC, AVIL, BAALC, BACH1, BIRC5, BRD2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 75 | `$find-subject-structural-intersecting-variants` | ALDH1A3, ALDH1L2, ALDH2, ALOX15B, ALOX5, ANK1, ANKRD26, APCDD1, APEX1, ART1, ASCL1, ASXL3, ATMIN, ATRIP, AURKC, AVIL, BAALC, BACH1, BIRC5, BRD2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 76 | `$find-subject-variants` | BRD9, CASR, CAV1, CCAR1, CCNB3, CD19, CDH2, CDH4, CDKN1C, CHTF8, CHUK, COL5A1, CTDNEP1, D2HGDH, DDX4, DDX41, DHX15, DIS3L2, DLL3, ECSIT (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 77 | `$find-subject-structural-intersecting-variants` | BRD9, CASR, CAV1, CCAR1, CCNB3, CD19, CDH2, CDH4, CDKN1C, CHTF8, CHUK, COL5A1, CTDNEP1, D2HGDH, DDX4, DDX41, DHX15, DIS3L2, DLL3, ECSIT (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 78 | `$find-subject-variants` | EGR2, EIF2B1, ELL2, ERCC1, ERCC6, ESCO1, ETV7, F3, FANCM, FGF1, FGF2, FGF5, FGF7, FGF8, FGF9, FOLH1, FOLR1, FOXN4, FSTL1, FZR1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 79 | `$find-subject-structural-intersecting-variants` | EGR2, EIF2B1, ELL2, ERCC1, ERCC6, ESCO1, ETV7, F3, FANCM, FGF1, FGF2, FGF5, FGF7, FGF8, FGF9, FOLH1, FOLR1, FOXN4, FSTL1, FZR1 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 80 | `$find-subject-variants` | GRB7, GSTO1, GSTP1, HDAC2, HFE, HSD17B2, HTATIP2, ID1, ING1, INTS6, IQGAP1, ITPKB, KAT7, KDM5D, KIAA1549, KLF2, KLF3, LGR5, LRP5, LRP6 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 81 | `$find-subject-structural-intersecting-variants` | GRB7, GSTO1, GSTP1, HDAC2, HFE, HSD17B2, HTATIP2, ID1, ING1, INTS6, IQGAP1, ITPKB, KAT7, KDM5D, KIAA1549, KLF2, KLF3, LGR5, LRP5, LRP6 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 82 | `$find-subject-variants` | LZTS1, MAGI2, MAL2, MAP3K21, MAP4K4, MASTL, MBD4, MBD6, MEF2D, MIDEAS, MLH3, MTHFD2, MYBL1, NQO1, NT5E, OGA, OGT, ONECUT2, PHF19 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 83 | `$find-subject-structural-intersecting-variants` | LZTS1, MAGI2, MAL2, MAP3K21, MAP4K4, MASTL, MBD4, MBD6, MEF2D, MIDEAS, MLH3, MTHFD2, MYBL1, NQO1, NT5E, OGA, OGT, ONECUT2, PHF19 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 84 | `$find-subject-variants` | PHLPP1, PHLPP2, POLH, POU2F2, POU3F2, POU3F4, PRKCB, PRMT5, PRPF8, PSMB2, PTPN14, RAD17, RELA, REV3L, RIOK2, RNASEH2A, RNASEH2B, RPS15, RPS6KB1, SF3B2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 85 | `$find-subject-structural-intersecting-variants` | PHLPP1, PHLPP2, POLH, POU2F2, POU3F2, POU3F4, PRKCB, PRMT5, PRPF8, PSMB2, PTPN14, RAD17, RELA, REV3L, RIOK2, RNASEH2A, RNASEH2B, RPS15, RPS6KB1, SF3B2 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 86 | `$find-subject-variants` | SFRP1, SFRP2, SLIT3, SQSTM1, SRP72, TACSTD2, TCL1B, TIGAR, TNFSF13, TOP2A, TRIB3, TYMS, UBE2A, UBTF, UCHL1, USP1, VEGFB, WEE1, WNT10A, WWOX (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 87 | `$find-subject-structural-intersecting-variants` | SFRP1, SFRP2, SLIT3, SQSTM1, SRP72, TACSTD2, TCL1B, TIGAR, TNFSF13, TOP2A, TRIB3, TYMS, UBE2A, UBTF, UCHL1, USP1, VEGFB, WEE1, WNT10A, WWOX (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 88 | `$find-subject-variants` | WWP1, XRCC1, XRCC3, YY1, ZBTB20, ZFP36L1, ZFP36L2, ZNF292, ZNF750 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |
| 89 | `$find-subject-structural-intersecting-variants` | WWP1, XRCC1, XRCC3, YY1, ZBTB20, ZFP36L1, ZFP36L2, ZNF292, ZNF750 (GRCh37 gene ranges) | Cancer-driver sweep: OncoKB Cancer Gene List genes typed ONCOGENE, TSG or ONCOGENE_AND_TSG | 0 variant(s) |

Round 1 variants (all somatic; VAF recorded as 0 for every variant, a placeholder):

| Variant (SPDI, GRCh37) | Gene | Change | On report? |
|---|---|---|---|
| NC_000001.10:27105549:C:T | ARID1A | NM_006015.6:c.5161C>T p.Arg1721* (stop-gained) | Yes |
| NC_000001.10:27106893::T | ARID1A | NM_006015.6:c.6504_6505insT p.Val2169Cysfs*56 (frameshift) | Yes |
| NC_000001.10:65301782:C:A | JAK1 | NM_002227.4/NM_001320923.2:c.3256G>T p.Ala1086Ser (missense, splice region) | Yes |
| NC_000016.9:347071:GG:GGG | AXIN1 | NM_003502.4:c.1939dup p.Arg647Profs*17 (frameshift) | No |
| NC_000007.13:41729606:G:A | INHBA | NM_002192.4:c.922C>T p.Arg308Cys (missense) | No |
| NC_000004.11:185320191:C:A | IRF2 | NM_002199.4:c.571G>T p.Glu191* (stop-gained) | No |

No structural variants (none possible: the test holds no CNV/SV types). Only 6 variants in 869 genes, so every variant found is listed in the synthesis; none were left to this log.

### Round 2: questions raised by round 1 and the evidence

| # | Operation | Query | Why | Result |
|---|---|---|---|---|
| 90 | `$find-study-metadata` | 20 ranges (minimum genes) | Studied/uncallable regions for absence claims | HTTP 503 (30 s host timeout); backed off |
| 91 | `$find-study-metadata` | 10 ranges | Retry with smaller batch | HTTP 503 ×2 (job stopped); switched to 3 ranges per request |
| 92 | `$find-study-metadata` | ARID1A range | Studied regions | 20 exon-level segments; c.5161 and c.6504 inside 27105513-27108595; no uncallable regions |
| 93–106 | `$find-study-metadata` | 3 ranges each: ARID2/AXIN1/JAK1; CTNNB1/TP53/IRF2; CD44/CDK9/CEACAM1; DNAJB1/EGFR/FGF19; FGF3/FGF4/FGFR1; FGFR2/FGFR4/GALNT14; KIF23/KRAS/NEDD9; PBLD/PRKACA/PROM1; PTPN6/RAF1/VEGFC; NTRK1/NTRK2/NTRK3; BRAF/RET/MLH1; MSH2/MSH6/PMS2; ERBB2/FGFR3/TERT promoter (NC_000005.9:1295100-1295300) | Minimum genes, genes stated unaltered (CTNNB1, TP53, BRAF), new findings (AXIN1, IRF2), and the TERT promoter (TERT mutated in 44% of HCC, Campani 2026) | 590 studied segments in all, all exon-level; uncallable regions: none. Covered: BRAF V600E (140453074-140453193), CTNNB1 exon 3 hotspots D32/S45 (41266016-41266244), TP53 R249 (7577498-7577608), AXIN1 R647, IRF2 E191. TERT promoter: only 1295100-1295183 studied; hotspots 1295228 (C228T) and 1295250 (C250T) NOT studied → not tested |
| 107 | `$find-subject-variants` | RYBP NC_000003.11:72420976-72496069; NHERF1 NC_000017.10:72744791-72765492; TRA NC_000014.8:22090056-23021074; TRB NC_000007.13:141998850-142510971; TRD NC_000014.8:22891536-22935568; TRG NC_000007.13:38279624-38407655; TERT promoter NC_000005.9:1295100-1295300 | The six OncoKB genes without server coordinates (fallback coordinates from HGNC → Ensembl/NCBI Gene, GRCh37; RYBP and SLC9A3R1/NHERF1 confirmed via Ensembl GRCh37 REST lookup today); TERT promoter question | 0 variants |
| 108 | `$find-subject-structural-intersecting-variants` | Same 7 ranges | Same | 0 variants |
| 109 | `$find-subject-molecular-consequences` | All 6 variants | Consequences from the server | Agree with VEP: ARID1A stop_gained (HIGH), ARID1A frameshift (HIGH, loss of function), AXIN1 c.1939dupC frameshift (HIGH, loss of function), IRF2 p.Glu191* stop_gained (HIGH, loss of function), JAK1 missense + splice region (MODERATE), INHBA p.Arg308Cys missense (MODERATE) |
| 110 | `$find-subject-tx-implications` | All 6 variants | Therapeutic implications | None returned (the report's own ruxolitinib/JAK1 level D implication is in the bundle only) |
| 111 | `$find-subject-dx-implications` | All 6 variants | Diagnostic implications / ClinVar | None returned |

Genome-server requests: 111 logged above, plus 869 coordinate lookups. Two rounds: round 1 (minimum + cancer-gene sweep), round 2 (studied regions, TERT promoter, unresolved genes, consequences and implications for the six variants). Round 2 raised no further genomic question: no trial or label names a finding not yet queried (CTNNB1 for NCT06811116 checked and absent — that trial also requires prior checkpoint inhibitor unless ineligible; NCT06915753 needs FGF19 IHC, not a genomic finding; fusions not testable; GPC3 and ARID1A protein need IHC).

### Driver assessment

| Variant | Gene role (OncoKB) | ClinVar / Hotspot | Driver? |
|---|---|---|---|
| ARID1A p.Arg1721* | TSG | Cancer Hotspots R1721 (13 tumors); ClinVar somatic Tier II | Yes (hotspot; truncating in TSG) |
| ARID1A p.Val2169fs | TSG | — | Yes (truncating in TSG) |
| AXIN1 p.Arg647fs | TSG | — | Yes (truncating in TSG) |
| IRF2 p.Glu191* | TSG | — | Yes (truncating in TSG) |
| JAK1 p.Ala1086Ser | ONCOGENE_AND_TSG | Not hotspot; not in ClinVar | No (missense, unknown effect; CIViC functional D: not gain of function) |
| INHBA p.Arg308Cys | ONCOGENE_AND_TSG | Not hotspot; ClinVar has only p.Arg308Ser (VUS), not this change | No (missense, unknown effect) |

Pathways (TCGA pathway files): AXIN1 in WNT; ARID1A, IRF2 not in any of the 10 pathways (labeled chromatin remodeling and interferon signaling in the synthesis).

### Reconciliation with the genetic report

- Reported ARID1A ×2 and JAK1 variants all found, same coordinates and HGVS.
- AXIN1 frameshift, IRF2 stop-gained and INHBA missense are in the genomic data but not on the report.
- Report specimen date 2026-07-22 vs server test date 2020-01-20.
- Report MSI-stable and TMB 1.3 can't be checked against the genomic data.
- The report's JAK1/ruxolitinib implication (level D, does not support) agrees with CIViC 10058.

## Knowledge-base, label, trial and literature lookups

All on 2026-11-16 unless stated.

### FDA labels (openFDA drug/label API)

| Query | Why | Result |
|---|---|---|
| `indications_and_usage:"hepatocellular"` (29 labels; Keytruda Qlex indication re-read after the always-search list was updated) | Standard options for HCC after sorafenib; therapies that might seem to apply | HCC indications extracted: Keytruda 2026-10-16 and Keytruda Qlex 2026-10-16 (HCC secondary to hepatitis B after prior systemic therapy other than PD-1/PD-L1); Stivarga 2026-02-09 and Cabometyx 2025-10-19 (previously treated with sorafenib); Cyramza 2026-06-09 (AFP ≥400 ng/mL after sorafenib); Tecentriq 2026-10-13 and Avastin 2025-01-06 (no prior systemic therapy); Opdivo 2026-08-12 and Yervoy 2026-08-20 (nivolumab + ipilimumab first-line, AND for HCC previously treated with sorafenib; first extraction missed the second indication, corrected after an independent fact-check); Opdivo Qvantig 2026-05-12 (monotherapy after sorafenib following IV nivolumab + ipilimumab); Imfinzi and Imjudo 2026-08-31 (unresectable HCC; HIMALAYA enrolled no prior systemic therapy); Lenvima 2026-10-14 (first-line); Nexavar 2023-08-28 (unresectable HCC) |
| Clinical studies sections of Keytruda, Stivarga, Cabometyx, Cyramza, Imfinzi | Ranking: trial design and survival results | KEYNOTE-394 OS 13.9 vs 13.0 mo, HR 0.78, ORR 11% vs 1.6%, DoR 23.9 mo, HBV viral load <2000 IU/mL required; RESORCE OS 10.6 vs 7.8 mo, HR 0.63, sorafenib tolerance ≥400 mg/day required; CELESTIAL OS 10.2 vs 8.0 mo, HR 0.76; REACH-2 OS 8.5 vs 7.3 mo, HR 0.71; HIMALAYA: no prior systemic therapy |
| Opdivo clinical studies (CheckMate-040 cohort 4), re-checked after fact-check | Nivolumab + ipilimumab after sorafenib | Single-arm, n=49, progressed on or intolerant to sorafenib, Child-Pugh A, 57% active HBV; ORR 33% (95% CI 20–48) by BICR RECIST, 4 CR; 56% of responses ≥12 months |
| `openfda.generic_name:"tazemetostat"` | Patient's current drug | Tazverik label 2026-04-29: epithelioid sarcoma and follicular lymphoma only; warning of secondary malignancies (1.7% MDS/AML/B-ALL at approval) |

### FDA and sponsor announcements (web)

| Query | Why | Result |
|---|---|---|
| "Ipsen tazemetostat discontinued clinical trials benefit risk 2026" | NCT05023655 terminated with this reason | Ipsen 2026-03-09: Tazverik withdrawn in all indications in all Ipsen markets; all tazemetostat trials and expanded access discontinued |
| "FDA Tazverik tazemetostat withdrawal Federal Register 2026" | FDA status | FDA alert, May 2026 (bulletin distributed 2026-05-12; FDA page current as of 2026-05-11): 18/318 (5.7%) hematologic second malignancies in SYMPHONY-1 vs 0 control; 3 deaths; cases after stopping; risks outweigh benefits; IDMC recommended stopping. Federal Register 2026-12367 (2026-06-22): NDA 211723 approval withdrawn |
| "FDA approves hepatocellular carcinoma 2026" | Approvals not yet in labels | No 2026 HCC approval found; rivoceranib + camrelizumab (first-line) received a complete response letter in July 2026 and was resubmitted |

### CIViC (GraphQL API, evidence items by molecular profile name)

| Query | Why | Result |
|---|---|---|
| ARID1A | Reported gene | 14 items; none rejected. Predictive: level B dasatinib does not support (ovarian clear cell); level D items (dasatinib, GSK126, PD-1/PD-L1, sorafenib in liver cancer 7328) — none level A/B for HCC therapy. Prognostic/oncogenic/functional items in other diseases |
| JAK1 | Reported gene | 16 items. A1086S: 10055 (functional, D, does not support gain of function), 10058 (predictive, D, HCC, ruxolitinib, does not support sensitivity) |
| AXIN1 | New finding | 0 items |
| IRF2 | New finding | 0 items |
| INHBA | New finding | 0 items |
| CIViC co-occurrence profiles (retrieved 2026-10-02) | Combinations | No profile includes ARID1A, AXIN1, IRF2, JAK1 or INHBA (file grep) |

### ClinVar (E-utilities, GRCh37 position search)

| Query | Result |
|---|---|
| 1[chr] AND 27105550[chrpos37] | 4530199 ARID1A c.5161C>T p.Arg1721Ter: somatic clinical impact Tier II – potential; no germline or oncogenicity classification |
| 1[chr] AND 27106890:27106896[chrpos37] | Only 4698986 ARID1A c.6507A>C (synonymous, likely benign); c.6504_6505insT not in ClinVar |
| 1[chr] AND 65301783[chrpos37] | No records (JAK1 c.3256G>T) |
| 4[chr] AND 185320192[chrpos37] | No records (IRF2 c.571G>T) |
| 7[chr] AND 41729607[chrpos37] | Only 4678850 INHBA c.922C>A p.Arg308Ser (VUS); c.922C>T not in ClinVar |
| 16[chr] AND 347060:347080[chrpos37] | Three nearby AXIN1 missense records (VUS, benign, unclassified); c.1939dup not in ClinVar |

### Cancer Hotspots (cancerhotspots.org API, single-residue hotspots by gene)

| Gene | Result |
|---|---|
| ARID1A | R1721 is a hotspot (13 tumors, q=1.5e-6) |
| JAK1 | Only L910; A1086 not a hotspot |
| AXIN1 | Only X652 (splice, liver); R647 not a hotspot |
| IRF2 | No hotspots |
| INHBA | No hotspots |

### OncoKB Cancer Gene List (downloaded 2026-11-06)

ARID1A TSG; AXIN1 TSG; IRF2 TSG; JAK1 ONCOGENE_AND_TSG; INHBA ONCOGENE_AND_TSG.

### Ensembl VEP (grch37.rest.ensembl.org, RefSeq transcripts)

| Variant | Result |
|---|---|
| 1:27105550 C>T | ARID1A NM_006015.6:c.5161C>T p.Arg1721Ter, stop_gained |
| 1:27106893 G>GT | ARID1A NM_006015.6:c.6504_6505insT p.Val2169CysfsTer56, frameshift |
| 1:65301783 C>A | JAK1 NM_002227.4 / NM_001320923.2 c.3256G>T p.Ala1086Ser, missense + splice region |
| 16:347071 C>CG | AXIN1 NM_003502.4:c.1939dup p.Arg647ProfsTer17, frameshift |
| 7:41729607 G>A | INHBA NM_002192.4:c.922C>T p.Arg308Cys, missense (SIFT deleterious, PolyPhen possibly damaging) |
| 4:185320192 C>A | IRF2 NM_002199.4:c.571G>T p.Glu191Ter, stop_gained |

### ClinicalTrials.gov (API v2)

| Query | Why | Result |
|---|---|---|
| Study NCT05023655 | Trial in the record | TERMINATED; primary completion 2026-03-06; why stopped: Ipsen discontinued all tazemetostat trials (unfavorable benefit-risk); single site Prisma Health, Greenville SC |
| Recruiting, `query.term=ARID1A` (8 studies) | ARID1A finding | NCT04104776 (tulmimetostat; excludes prior EZH2 inhibitor and active HBV/HCV) — not eligible; NCT05490472 (JAB-2485; ARID1A-mutant expansion; excludes active HBV) — listed with caveat; NCT06811116 (needs β-catenin mutation) — not eligible; others gastroesophageal/gyn or no US site |
| Recruiting interventional, `query.cond=hepatocellular carcinoma`, US (121 studies) | Disease setting | Eligibility reviewed for advanced-HCC treatment trials after prior systemic therapy: listed NCT05440708, NCT05797805, NCT06600321, NCT06858813, NCT07042919, NCT07174570, NCT04864054, NCT06764316, NCT06915753 (re-read after fact-check: Part A HCC after standard care with no genomic requirement; Part B HCC cohort needs archival tissue for central FGF19 IHC). Not eligible: NCT06726161 (prior interventional study participation), NCT05733598 / NCT07219459 / NCT05168163 (prior PD-(L)1 required), NCT06294548 (no prior systemic therapy), NCT07175441 (first-line STRIDE), NCT07166406 (macrovascular invasion required) |
| Recruiting, `query.term=IRF2` | New IRF2 finding | 0 studies |
| Recruiting, `query.term=AXIN1` | New AXIN1 finding | 3 studies, none with a US site or relevant (RK-582 colorectal; non-oncology; renal) |
| Tumor-agnostic trial survey (2026-11-16) | Checklist | ARID1A → NCT04104776 (not eligible, prior EZH2 inhibitor). JAK1, AXIN1, IRF2, INHBA not in survey. No other match. NCT05490472 (JAB-2485, ARID1A-mutant solid tumors) found here is not in the survey: reported to the coordinator |

### PubMed (E-utilities, publication date 2025-10-19 to 2026-11-16)

| Search | Count | Used |
|---|---|---|
| ARID1A AND (hepatocellular OR liver cancer) | 25 | Campani 2026 J Hepatol (French genomic medicine initiative) |
| tazemetostat AND ARID1A | 1 | Computational study only |
| (EZH2 inhibitor OR tazemetostat) AND solid tumor AND clinical trial | 3 | Gounder 2026 Nat Commun (phase 2 basket) |
| ARID1A AND (ATR inhibitor OR ceralasertib OR camonsertib) | 9 | Zhu 2026 Clin Cancer Res (ceralasertib phase 2) |
| hepatocellular carcinoma AND (second-line OR after sorafenib) AND (phase 3 OR randomized) | 34 | KEYNOTE-224 long-term follow-up (Clin Cancer Res 2026) reviewed; consistent with labels; network meta-analyses not used |
| JAK1 AND hepatocellular carcinoma | 15 | Nothing changing the JAK1 reading |
| IRF2 AND hepatocellular carcinoma | 1 | Genomic subtype study; no therapeutic evidence |
| EZH2 inhibitor AND hepatocellular carcinoma | 11 | Preclinical/reviews only |
| AXIN1 AND hepatocellular carcinoma | 12 | iScience 2026 (AXIN1 mutations drive lower WNT activity) |
| (beta-catenin OR CTNNB1 OR WNT) AND hepatocellular carcinoma AND (immunotherapy OR immune checkpoint) AND (response OR resistance) | 52 | CheckMate 459 biomarker analysis (JHEP Rep 2026) |

Search window: 2025-10-19 to 2026-11-16. The window starts at the earliest of: the earliest label cited (Cabometyx, 2025-10-19), the co-occurrence knowledge base retrieval (2026-10-02) and the AI training cutoff (June 2026).
