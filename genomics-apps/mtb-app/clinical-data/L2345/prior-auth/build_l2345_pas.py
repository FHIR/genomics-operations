"""Build a synthetic Da Vinci PAS (STU 2.0.1) prior authorization request and response for L2345.

The request asks for datopotamab deruxtecan (Datroway), one of the FDA-labeled options in L2345's case
synthesis. Clinical justification is copied from L2345-transaction-bundle.json (the patient's mCODE record)
plus his EGFR L858R variant from FHIR Genomics Operations. Coverage, payer, NPIs, member id and the payer's
decision are invented; see README.md.

Run from any directory:  python3 build_l2345_pas.py
"""

import copy
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_BUNDLE = HERE.parent / 'L2345-transaction-bundle.json'

PAS = 'http://hl7.org/fhir/us/davinci-pas/StructureDefinition/'
X12 = 'https://codesystem.x12.org/005010/'
HCPCS = 'http://www.cms.gov/Medicare/Coding/HCPCSReleaseCodeSets'
ICD10CM = 'http://hl7.org/fhir/sid/icd-10-cm'
NPI = 'http://hl7.org/fhir/sid/us-npi'
LOINC = 'http://loinc.org'
FULL_URL = 'http://elimu.io/fhir/'
DEMO_TAG = {'system': 'http://elimu.io/fhir/tags', 'code': 'mtb-demo', 'display': 'MTB app synthetic demo data'}

# Dates: request submitted today; six 21-day cycles starting a week later
CREATED = '2026-11-02T10:15:00-05:00'
RESPONDED = '2026-11-04T14:40:00-05:00'
SERVICE_PERIOD = {'start': '2026-11-09', 'end': '2027-02-22'}
AUTH_PERIOD = {'start': '2026-11-09', 'end': '2027-03-08'}

# Dose per label: 6 mg/kg IV every 3 weeks, max 540 mg; weight 74.2 kg on 2026-10-08
WEIGHT_KG = 74.2
DOSE_MG = min(round(6 * WEIGHT_KG), 540)
CYCLES = 6
UNITS = DOSE_MG * CYCLES  # J9011 is billed per 1 mg

# Clinical resources from the patient's record that justify the request
SUPPORTING_IDS = [
    'mtb-l2345-primary-cancer',
    'mtb-l2345-liver-metastases',
    'mtb-l2345-tnm-clinical-stage-group',
    'mtb-l2345-disease-status-2026-09-28',
    'mtb-l2345-ct-cap-2026-09-28',
    'mtb-l2345-path-liver-biopsy',
    'mtb-l2345-medrx-osimertinib',
    'mtb-l2345-medrx-carboplatin',
    'mtb-l2345-medrx-pemetrexed',
    'mtb-l2345-molecular-2023-06-09',
    'mtb-l2345-ecog-2026-10-08',
    'mtb-l2345-weight-2026-10-08',
]


def npi(base9: str) -> str:
    """A 10-digit NPI with a valid Luhn check digit (computed over the 80840 prefix)."""
    digits = [int(d) for d in '80840' + base9]
    total = 0
    for index, digit in enumerate(reversed(digits)):
        if index % 2 == 0:
            digit *= 2
            if digit > 9:
                digit -= 9
        total += digit
    return base9 + str((10 - total % 10) % 10)


def ref(resource):
    return {'reference': f"{resource['resourceType']}/{resource['id']}"}


def entry(resource):
    return {'fullUrl': f"{FULL_URL}{resource['resourceType']}/{resource['id']}", 'resource': resource}


def with_profile(resource, profile):
    resource = copy.deepcopy(resource)
    meta = resource.setdefault('meta', {})
    meta['profile'] = [PAS + profile] + [p for p in meta.get('profile', []) if p != PAS + profile]
    return resource


def x12(code_list, code, display=None):
    coding = {'system': X12 + code_list, 'code': code}
    if display:
        coding['display'] = display
    return {'coding': [coding]}


source = {e['resource']['id']: e['resource'] for e in json.loads(SOURCE_BUNDLE.read_text())['entry']}

# --- Administrative data (synthetic) ---------------------------------------------------------------

insurer = {
    'resourceType': 'Organization', 'id': 'mtb-l2345-insurer',
    'meta': {'profile': [PAS + 'profile-insurer'], 'tag': [DEMO_TAG]},
    'identifier': [{'system': 'http://elimu.io/fhir/synthetic-health-plan/payer-id', 'value': 'SHP01'}],
    'active': True,
    'type': [x12('98', 'PR', 'Payer')],
    'name': 'Synthetic Health Plan (demo)',
}

requestor = {
    'resourceType': 'Organization', 'id': 'mtb-l2345-requestor',
    'meta': {'profile': [PAS + 'profile-requestor'], 'tag': [DEMO_TAG]},
    'identifier': [{'system': NPI, 'value': npi('730012345')}],
    'active': True,
    'type': [x12('98', 'X3')],
    'name': 'MTB Demo Cancer Center',
    'address': [{'line': ['100 Demo Way'], 'city': 'Springfield', 'state': 'IL', 'postalCode': '62701', 'country': 'US'}],
}

practitioner = with_profile(source['mtb-l2345-oncologist'], 'profile-practitioner')
practitioner['identifier'] = [{'system': NPI, 'value': npi('520098765')}] + practitioner.get('identifier', [])

practitioner_role = {
    'resourceType': 'PractitionerRole', 'id': 'mtb-l2345-oncologist-role',
    'meta': {'profile': [PAS + 'profile-practitionerrole'], 'tag': [DEMO_TAG]},
    'practitioner': ref(practitioner),
    'organization': ref(requestor),
    'code': [{'coding': [{'system': 'http://nucc.org/provider-taxonomy', 'code': '207RX0202X',
                          'display': 'Medical Oncology Physician'}]}],
}

patient = with_profile(source['mtb-l2345'], 'profile-beneficiary')
patient['identifier'] = patient.get('identifier', []) + [{
    'type': {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/v2-0203', 'code': 'MB', 'display': 'Member Number'}]},
    'system': 'http://elimu.io/fhir/synthetic-health-plan/member-id',
    'value': 'SHP884120577',
}]

coverage = {
    'resourceType': 'Coverage', 'id': 'mtb-l2345-coverage',
    'meta': {'profile': [PAS + 'profile-coverage'], 'tag': [DEMO_TAG]},
    'status': 'active',
    'subscriberId': 'SHP884120577',
    'beneficiary': ref(patient),
    'relationship': {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/subscriber-relationship', 'code': 'self'}]},
    'period': {'start': '2026-01-01', 'end': '2026-12-31'},
    'payor': [ref(insurer)],
    'class': [
        {'type': {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/coverage-class', 'code': 'group'}]},
         'value': 'GRP-20410', 'name': 'Demo Employer Group'},
        {'type': {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/coverage-class', 'code': 'plan'}]},
         'value': 'PPO-GOLD', 'name': 'Gold PPO'},
    ],
}

# --- The requested drug -----------------------------------------------------------------------------

medication_request = {
    'resourceType': 'MedicationRequest', 'id': 'mtb-l2345-medrx-datopotamab',
    'meta': {'profile': [PAS + 'profile-medicationrequest'], 'tag': [DEMO_TAG]},
    'status': 'draft',
    'intent': 'order',
    'medicationCodeableConcept': {
        'coding': [{'system': 'http://hl7.org/fhir/sid/ndc', 'code': '65597-801-01'}],
        'text': 'datopotamab deruxtecan-dlnk (Datroway) 100 mg vial',
    },
    'subject': ref(patient),
    'authoredOn': CREATED[:10],
    'requester': ref(practitioner),
    'reasonReference': [{'reference': 'Condition/mtb-l2345-primary-cancer'}],
    'dosageInstruction': [{
        'text': f'6 mg/kg ({DOSE_MG} mg at {WEIGHT_KG} kg) IV infusion every 3 weeks, {CYCLES} cycles, then reassess',
        'timing': {'repeat': {'count': CYCLES, 'frequency': 1, 'period': 3, 'periodUnit': 'wk'}},
        'route': {'coding': [{'system': 'http://snomed.info/sct', 'code': '47625008', 'display': 'Intravenous route'}]},
        'doseAndRate': [{'doseQuantity': {'value': DOSE_MG, 'unit': 'mg', 'system': 'http://unitsofmeasure.org', 'code': 'mg'}}],
    }],
}

# --- Genomic evidence from FHIR Genomics Operations ---------------------------------------------------
# $find-subject-variants (subject L2345, EGFR range) returned Observation dv-69851b51cd18150898d85c1d with SPDI
# NC_000007.14:55191821:T:G; $find-subject-molecular-consequences gave its c. and p. HGVS. Queried 2026-11-02.

egfr_variant = {
    'resourceType': 'Observation', 'id': 'mtb-l2345-egfr-l858r',
    'meta': {'tag': [DEMO_TAG]},
    'status': 'final',
    'category': [{'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/observation-category', 'code': 'laboratory'}]}],
    'code': {'coding': [{'system': LOINC, 'code': '69548-6', 'display': 'Genetic variant assessment'}]},
    'subject': ref(patient),
    'valueCodeableConcept': {'coding': [{'system': LOINC, 'code': 'LA9633-4', 'display': 'Present'}]},
    'note': [{'text': 'Source: FHIR Genomics Operations (fhir-gen-ops.herokuapp.com), variant dv-69851b51cd18150898d85c1d.'}],
    'component': [
        {'code': {'coding': [{'system': LOINC, 'code': '48018-6', 'display': 'Gene studied [ID]'}]},
         'valueCodeableConcept': {'coding': [{'system': 'http://www.genenames.org/geneId', 'code': 'HGNC:3236', 'display': 'EGFR'}]}},
        {'code': {'coding': [{'system': LOINC, 'code': '81252-9', 'display': 'Discrete genetic variant'}]},
         'valueCodeableConcept': {'coding': [{'system': 'https://api.ncbi.nlm.nih.gov/variation/v0/', 'code': 'NC_000007.14:55191821:T:G'}]}},
        {'code': {'coding': [{'system': LOINC, 'code': '48004-6', 'display': 'DNA change (c.HGVS)'}]},
         'valueCodeableConcept': {'text': 'NM_005228.5:c.2573T>G'}},
        {'code': {'coding': [{'system': LOINC, 'code': '48005-3', 'display': 'Amino acid change (pHGVS)'}]},
         'valueCodeableConcept': {'text': 'NP_005219.2:p.Leu858Arg'}},
    ],
}

supporting = [copy.deepcopy(source[i]) for i in SUPPORTING_IDS] + [egfr_variant]

# --- The PAS Claim ----------------------------------------------------------------------------------

necessity = (
    f'69-year-old never-smoker with stage IVB EGFR L858R-mutated adenocarcinoma of the right lower lobe '
    f'(C34.31), diagnosed 2023-05. First line: osimertinib with carboplatin and pemetrexed from 2023-06, then '
    f'pemetrexed maintenance. Progression on 2026-09-28 with new liver metastases (biopsy-confirmed '
    f'adenocarcinoma, no small cell transformation) and rising CEA. ECOG 1. He has received both prior '
    f'EGFR-directed therapy and platinum-based chemotherapy, as required by the Datroway indication for EGFR-'
    f'mutated NSCLC. Requesting datopotamab deruxtecan 6 mg/kg ({DOSE_MG} mg) IV every 3 weeks for {CYCLES} '
    f'cycles, then reassessment.'
)

claim = {
    'resourceType': 'Claim', 'id': 'mtb-l2345-pas-claim',
    'meta': {'profile': [PAS + 'profile-claim'], 'tag': [DEMO_TAG]},
    'identifier': [{'system': 'http://elimu.io/fhir/mtb-demo/patient-event-trace-number', 'value': 'L2345-PA-0001'}],
    'status': 'active',
    'type': {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/claim-type', 'code': 'professional'}]},
    'use': 'preauthorization',
    'patient': ref(patient),
    'created': CREATED,
    'insurer': ref(insurer),
    'provider': ref(requestor),
    'priority': {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/processpriority', 'code': 'normal'}]},
    'careTeam': [{
        'extension': [{'url': PAS + 'extension-careTeamClaimScope', 'valueBoolean': True}],
        'sequence': 1,
        'provider': ref(practitioner_role),
    }],
    'supportingInfo': [
        {'sequence': 1, 'category': {'coding': [{'system': 'http://hl7.org/fhir/us/davinci-pas/CodeSystem/PASTempCodes', 'code': 'freeFormMessage'}]},
         'valueString': necessity},
    ] + [
        {'sequence': n, 'category': {'coding': [{'system': 'http://hl7.org/fhir/us/davinci-pas/CodeSystem/PASTempCodes', 'code': 'additionalInformation'}]},
         'valueReference': ref(resource)}
        for n, resource in enumerate(supporting, start=2)
    ],
    'diagnosis': [
        {'sequence': n, 'diagnosisCodeableConcept': {'coding': [{'system': ICD10CM, 'code': code, 'display': display}]}}
        for n, (code, display) in enumerate([
            ('C34.31', 'Malignant neoplasm of lower lobe, right bronchus or lung'),
            ('C78.7', 'Secondary malignant neoplasm of liver and intrahepatic bile duct'),
            ('C79.51', 'Secondary malignant neoplasm of bone'),
        ], start=1)
    ],
    'insurance': [{'sequence': 1, 'focal': True, 'coverage': ref(coverage)}],
    'item': [{
        'extension': [
            {'url': PAS + 'extension-itemTraceNumber',
             'valueIdentifier': {'system': 'http://elimu.io/fhir/mtb-demo/item-trace-number', 'value': 'L2345-PA-0001-1'}},
            {'url': PAS + 'extension-serviceItemRequestType', 'valueCodeableConcept': x12('1525', 'HS', 'Health Services Review')},
            {'url': PAS + 'extension-certificationType', 'valueCodeableConcept': x12('1322', 'I', 'Initial')},
            {'url': PAS + 'extension-requestedService', 'valueReference': ref(medication_request)},
        ],
        'sequence': 1,
        'careTeamSequence': [1],
        'diagnosisSequence': [1, 2, 3],
        'category': x12('1365', '78', 'Chemotherapy'),
        'productOrService': {'coding': [{'system': HCPCS, 'code': 'J9011',
                                         'display': 'Injection, datopotamab deruxtecan-dlnk, 1 mg'}]},
        'servicedPeriod': SERVICE_PERIOD,
        'locationCodeableConcept': {'coding': [{'system': 'https://www.cms.gov/Medicare/Coding/place-of-service-codes/Place_of_Service_Code_Set',
                                                'code': '11', 'display': 'Office'}]},
        'quantity': {'value': UNITS, 'unit': 'units (1 mg each)'},
    }],
}

request_bundle = {
    'resourceType': 'Bundle', 'id': 'mtb-l2345-pas-request',
    'meta': {'profile': ['http://hl7.org/fhir/us/davinci-pas/StructureDefinition/profile-pas-request-bundle'], 'tag': [DEMO_TAG]},
    'identifier': {'system': 'http://elimu.io/fhir/mtb-demo/submitter-transaction-identifier', 'value': 'L2345-PA-0001'},
    'type': 'collection',
    'timestamp': CREATED,
    'entry': [entry(r) for r in [claim, patient, coverage, insurer, requestor, practitioner_role, practitioner,
                                 medication_request] + supporting],
}

# --- The payer's response (synthetic decision) -----------------------------------------------------------

claim_response = {
    'resourceType': 'ClaimResponse', 'id': 'mtb-l2345-pas-claimresponse',
    'meta': {'profile': [PAS + 'profile-claimresponse'], 'tag': [DEMO_TAG]},
    'identifier': claim['identifier'],
    'status': 'active',
    'type': claim['type'],
    'use': 'preauthorization',
    'patient': ref(patient),
    'created': RESPONDED,
    'insurer': ref(insurer),
    'requestor': ref(requestor),
    'request': ref(claim),
    'outcome': 'complete',
    'disposition': 'Certified in total',
    'preAuthRef': 'SHP-PA-2026-110457',
    'preAuthPeriod': AUTH_PERIOD,
    'item': [{
        'extension': [
            {'url': PAS + 'extension-itemTraceNumber',
             'valueIdentifier': claim['item'][0]['extension'][0]['valueIdentifier']},
            {'url': PAS + 'extension-itemPreAuthIssueDate', 'valueDate': RESPONDED[:10]},
            {'url': PAS + 'extension-itemPreAuthPeriod', 'valuePeriod': AUTH_PERIOD},
            {'url': PAS + 'extension-itemAuthorizedDetail', 'extension': [
                {'url': 'productOrServiceCode', 'valueCodeableConcept': claim['item'][0]['productOrService']},
                {'url': 'quantity', 'valueQuantity': {'value': UNITS}},
            ]},
            {'url': PAS + 'extension-itemAuthorizedProvider', 'extension': [
                {'url': 'provider', 'valueReference': ref(practitioner)},
            ]},
        ],
        'itemSequence': 1,
        'adjudication': [{
            'extension': [{'url': PAS + 'extension-reviewAction', 'extension': [
                {'url': 'number', 'valueString': 'SHP-PA-2026-110457'},
                {'url': PAS + 'extension-reviewActionCode', 'valueCodeableConcept': x12('306', 'A1', 'Certified in total')},
            ]}],
            'category': {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/adjudication', 'code': 'submitted'}]},
        }],
    }],
    'processNote': [{
        'number': 1,
        'type': 'display',
        'text': (f'Approved: {UNITS} units of J9011 ({CYCLES} cycles of {DOSE_MG} mg), '
                 f'{AUTH_PERIOD["start"]} to {AUTH_PERIOD["end"]}. Criteria met: EGFR-mutated metastatic NSCLC; '
                 'progression after an EGFR TKI (osimertinib) and platinum-based chemotherapy; ECOG 0-2. '
                 'Continuation beyond this period requires a new request with imaging showing no disease progression.'),
    }],
}

response_bundle = {
    'resourceType': 'Bundle', 'id': 'mtb-l2345-pas-response',
    'meta': {'profile': ['http://hl7.org/fhir/us/davinci-pas/StructureDefinition/profile-pas-response-bundle'], 'tag': [DEMO_TAG]},
    'identifier': {'system': 'http://elimu.io/fhir/synthetic-health-plan/response-id', 'value': 'SHP-RESP-2026-110457'},
    'type': 'collection',
    'timestamp': RESPONDED,
    'entry': [entry(r) for r in [claim_response, patient, insurer, requestor, practitioner]],
}

for name, bundle in [('L2345-pas-request-bundle.json', request_bundle), ('L2345-pas-response-bundle.json', response_bundle)]:
    (HERE / name).write_text(json.dumps(bundle, indent=2) + '\n')
    print(f'wrote {name}: {len(bundle["entry"])} entries')
print(f'dose {DOSE_MG} mg x {CYCLES} cycles = {UNITS} units of J9011')
