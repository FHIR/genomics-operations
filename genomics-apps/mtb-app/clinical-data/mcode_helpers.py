"""Shared helpers for the build_<MRN>.py scripts that generate mCODE transaction bundles.

Each build script runs in its own process: resources passed to add() collect in the
module-level `entries` list, and write_bundle() writes them out as one transaction Bundle.
"""
import json
import re
from datetime import date

MCODE = 'http://hl7.org/fhir/us/mcode/StructureDefinition/'
USCORE = 'http://hl7.org/fhir/us/core/StructureDefinition/'
SCT = 'http://snomed.info/sct'
LOINC = 'http://loinc.org'
RXN = 'http://www.nlm.nih.gov/research/umls/rxnorm'
ICD10 = 'http://hl7.org/fhir/sid/icd-10-cm'
ICDO3 = 'http://terminology.hl7.org/CodeSystem/icd-o-3'
UCUM = 'http://unitsofmeasure.org'
TAG = {'system': 'http://elimu.io/fhir/tags', 'code': 'mtb-demo', 'display': 'MTB app synthetic demo data'}

OBS_CAT = 'http://terminology.hl7.org/CodeSystem/observation-category'
COND_CLINICAL = 'http://terminology.hl7.org/CodeSystem/condition-clinical'
COND_VER = 'http://terminology.hl7.org/CodeSystem/condition-ver-status'
COND_CAT = 'http://terminology.hl7.org/CodeSystem/condition-category'

LEFT = ('7771000', 'Left')
RIGHT = ('24028007', 'Right')

entries = []


def cc(system, code, display, text=None):
    c = {'coding': [{'system': system, 'code': code, 'display': display}]}
    if text:
        c['text'] = text
    return c


def add(resource, profiles=()):
    rid = resource['id']
    resource = {'resourceType': resource.pop('resourceType'), 'id': rid,
                'meta': {**({'profile': list(profiles)} if profiles else {}), 'tag': [TAG]}, **resource}
    entries.append({
        'fullUrl': f"https://hapi.fhir.org/baseR4/{resource['resourceType']}/{rid}",
        'resource': resource,
        'request': {'method': 'PUT', 'url': f"{resource['resourceType']}/{rid}"},
    })
    return {'reference': f"{resource['resourceType']}/{rid}"}


def ext(name, value_key, value):
    return {'url': MCODE + name, value_key: value}


def qty(value, unit, code=None):
    return {'value': value, 'unit': unit, 'system': UCUM, 'code': code or unit}


def obs_cat(code, display):
    return [cc(OBS_CAT, code, display)]


def condition_base(status='active'):
    return {
        'clinicalStatus': cc(COND_CLINICAL, status, status.capitalize()),
        'verificationStatus': cc(COND_VER, 'confirmed', 'Confirmed'),
        'category': [cc(COND_CAT, 'problem-list-item', 'Problem List Item')],
    }


def laterality(code, display):
    return [{'url': MCODE + 'mcode-laterality-qualifier', 'valueCodeableConcept': cc(SCT, code, display)}]


def specimen_type(v2_code, v2_display, sct_code, sct_display):
    return {'coding': [{'system': 'http://terminology.hl7.org/CodeSystem/v2-0487', 'code': v2_code, 'display': v2_display},
                       {'system': SCT, 'code': sct_code, 'display': sct_display}]}


def relink(node, id_map):
    """Rewrite every reference to a relative `Type/new-id`, using id_map {(type, old id): new id}.

    Handles absolute references (including ones with another server's base URL prefixed)."""
    if isinstance(node, dict):
        for k, v in node.items():
            if k == 'reference' and isinstance(v, str):
                rtype, rid = v.rstrip('/').split('/')[-2:]
                node[k] = f'{rtype}/{id_map[(rtype, rid)]}'
            else:
                relink(v, id_map)
    elif isinstance(node, list):
        for v in node:
            relink(v, id_map)


DATE_RE = re.compile(r'\b(19|20)\d\d-\d\d-\d\d\b')


def shift_dates(node, shift):
    """Move every date in node by the timedelta `shift`: date and dateTime values, dates inside
    ids, references and text, and year-only *DateTime values (by the nearest whole number of years)."""
    if isinstance(node, dict):
        for k, v in node.items():
            if isinstance(v, str):
                if k.endswith('DateTime') and re.fullmatch(r'(19|20)\d\d', v):
                    node[k] = str(int(v) + round(shift.days / 365.25))
                else:
                    node[k] = DATE_RE.sub(lambda m: (date.fromisoformat(m.group()) + shift).isoformat(), v)
            else:
                shift_dates(v, shift)
    elif isinstance(node, list):
        for v in node:
            shift_dates(v, shift)


def write_bundle(out, bundle_id):
    bundle = {'resourceType': 'Bundle', 'id': bundle_id, 'meta': {'tag': [TAG]},
              'type': 'transaction', 'entry': entries}
    with open(out, 'w', newline='\n') as f:
        json.dump(bundle, f, indent=2, ensure_ascii=False)
        f.write('\n')
    print(f'{len(entries)} entries written to {out}')
