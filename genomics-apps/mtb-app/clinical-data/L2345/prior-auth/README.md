# L2345 prior authorization example (Da Vinci PAS)

A synthetic prior authorization request and payer response for Larry Lung (MRN L2345), built to show how the MTB app's clinical and genomic data could justify a new treatment. Not for clinical or billing use.

**Files**
- [build_l2345_pas.py](build_l2345_pas.py): generator (edit this, not the JSON).
- [L2345-pas-request-bundle.json](L2345-pas-request-bundle.json): PAS request bundle (Claim with `use` = preauthorization, plus everything it references).
- [L2345-pas-response-bundle.json](L2345-pas-response-bundle.json): PAS response bundle (ClaimResponse).

Profiles: Da Vinci PAS STU 2.0.1 (`hl7.fhir.us.davinci-pas#2.0.1`), FHIR R4.

## The request

**Datopotamab deruxtecan (Datroway)**, one of the FDA-labeled options in L2345's case synthesis (`public/data/synthesis/L2345.md`). Its NSCLC indication requires EGFR-mutated disease and prior EGFR-directed therapy and platinum-based chemotherapy, all of which he has. It also avoids repeating his prior chemotherapy backbone (the caveat on amivantamab + carboplatin + pemetrexed).

- **Dose:** 6 mg/kg IV every 3 weeks (label maximum 540 mg). At 74.2 kg (2026-10-08) that is 445 mg per cycle.
- **Quantity:** 6 cycles, 2,670 units of HCPCS **J9011** (billed per 1 mg), from 2026-11-09.
- **Diagnoses:** ICD-10-CM C34.31 (primary), C78.7 (liver metastases), C79.51 (bone metastases).
- **Item category:** X12 service type 78, Chemotherapy. Place of service 11, Office.

## Where each part comes from

| Part | Source |
|---|---|
| Patient, primary cancer, liver metastases, stage IVB, progression 2026-09-28, CT and liver biopsy reports, prior osimertinib, carboplatin and pemetrexed, 2023 molecular report, ECOG 1, weight | **Real demo data:** copied unchanged from `../L2345-transaction-bundle.json` (the patient's mCODE record on HAPI) |
| EGFR L858R variant (SPDI `NC_000007.14:55191821:T:G`, `NM_005228.5:c.2573T>G`, `NP_005219.2:p.Leu858Arg`) | **Real demo data:** FHIR Genomics Operations, `$find-subject-variants` and `$find-subject-molecular-consequences`, queried 2026-11-02 |
| Free-text medical necessity summary | Written from the two sources above |
| Payer (Synthetic Health Plan), coverage, member id SHP884120577, group and plan | **Invented** |
| Requesting organization (MTB Demo Cancer Center), NPIs (Luhn-valid, fictional), the oncologist's PractitionerRole | **Invented** |
| NDC 65597-0801-01, J9011, dosing | From public billing guidance and the label summary (verify before reuse) |
| The decision: certified in total (X12 306 A1), authorization SHP-PA-2026-110457, valid 2026-11-09 to 2027-03-08, renewal requires imaging without progression | **Invented** (a plausible approval: the request meets the label criteria) |

The request is dated 2026-11-02 and the response 2026-11-04, within the 7-day standard decision time required by the CMS Interoperability and Prior Authorization Final Rule (CMS-0057-F).

## Not modeled
- CRD (asking the payer whether prior authorization is needed) and DTR (the payer's questionnaire). The justification here is free text plus referenced resources, not a QuestionnaireResponse.
- Submission (`Claim/$submit`), pended responses, status inquiry and subscriptions.
