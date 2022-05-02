from collections import OrderedDict
from uuid import uuid4
import pyliftover
import requests
from datetime import datetime
import pymongo
from flask import abort
from itertools import groupby
import os
import re

# Environment Variables
MONGODB_USERNAME = os.environ.get('MONGODB_USERNAME')
MONGODB_PASSWORD = os.environ.get('MONGODB_PASSWORD')

# MongoDB Clients
client = pymongo.MongoClient(f"mongodb+srv://{MONGODB_USERNAME}:{MONGODB_PASSWORD}@cluster0.pe0bv.mongodb.net/FHIRGenomicsData")
utilities_client = pymongo.MongoClient(f"mongodb+srv://{MONGODB_USERNAME}:{MONGODB_PASSWORD}@cluster0.8ianr.mongodb.net/UtilitiesData")

# Databases
db = client.FHIRGenomicsData
utilities_db = utilities_client.UtilitiesData

# Collections
patients_db = db.Patients
variants_db = db.Variants
tests_db = db.Tests
genotypes_db = db.Genotypes
clinvar_db = db.Clinvar
civic_db = db.CIViC
pharmGKB_db = db.PharmGKB
beds = db.BEDs

chromosomes_data = utilities_db.Chromosomes
genes_data = utilities_db.Genes
transcripts_data = utilities_db.Transcripts
exons_data = utilities_db.Exons
proteins_data = utilities_db.Proteins

# Constants
SVTYPE_TO_DNA_CHANGE_TYPE = {
    'CNV': ['SO:0001019', 'copy_number_variation'],
    'DUP': ['SO:1000035', 'duplication'],
    'INV': ['SO:1000036', 'inversion'],
    'DEL': ['SO:0000159', 'deletion'],
    'INS': ['SO:0000667', 'insertion']
}

DNA_CHANGE_TYPE_TO_CODE = {
    'SNV': 'SO:0001483',
    'MNV': 'SO:0002007',
    'INS': 'SO:0000667',
    'DEL': 'SO:0000159',
    'INV': 'SO:1000036',
    'DUP': 'SO:0001019',
    'CNV': 'SO:0001019',
    'delins': 'SO:1000032',
    'copy_number_variation': 'SO:0001019'
}

GENOMIC_SOURCE_CLASS_TO_CODE = {
    'GERMLINE': 'LA6683-2',
    'SOMATIC': 'LA6684-0'
}

GENOMIC_BUILD_TO_CODE = {
    'GRCh37': 'LA14029-5',
    'GRCh38': ' LA26806-2'
}

ALLELIC_STATE_TO_CODE = {
    'HETEROZYGOUS': 'LA6706-1',
    'HOMOZYGOUS': 'LA6705-3',
    'HEMIZYGOUS': 'LA6707-9'
}

CLIN_SIG_TO_CODE = {
    'Pathogenic': 'LA6668-3',
    'Likely pathogenic': 'LA26332-9',
    'Uncertain significance': 'LA26333-7',
    'Likely benign': 'LA26334-5',
    'Benign': 'LA6675-8'
}

OPERATOR_CODE_TO_OPERATOR = {
    'eq': '$eq',
    'ne': '$ne',
    'lt': '$lt',
    'gt': '$gt',
    'ge': '$gte',
    'le': '$lte',
    'sa': '$gte',
    'eb': '$lte'
}

SUPPORTED_GENE_SYSTEM_URLS = r'^https?:\/\/www\.genenames.org\/geneId$'

SUPPORTED_DATE_FORMAT = '%Y-%m-%d'

SUPPORTED_GENOMIC_SOURCE_CLASSES = ['germline', 'somatic']

NCBI_VARIATION_SERVICES_BASE_URL = 'https://api.ncbi.nlm.nih.gov/variation/v0/'

CHROMOSOME_CSV_FILE = 'app/_Dict_Chromosome.csv'

# Utility Functions
def key_func(k):
    return k['CHROMOSOME']['RefSeq']


def is_overlaping(a, b):
    if b['RANGE']['L'] >= a['RANGE']['L'] and b['RANGE']['L'] <= a['RANGE']['H']:
        return True
    else:
        return False


def merge(arr):
    arr.sort(key=lambda x: x['RANGE']['L'])
    merged_list = []
    merged_list.append(arr[0])

    for i in range(1, len(arr)):
        pop_element = merged_list.pop()

        if is_overlaping(pop_element, arr[i]):
            pop_element['RANGE']['H'] = max(pop_element['RANGE']['H'], arr[i]['RANGE']['H'])
            merged_list.append(pop_element)
        else:
            merged_list.append(pop_element)
            merged_list.append(arr[i])

    return merged_list


def merge_ranges(ranges):  
    ranges = sorted(ranges, key=key_func)

    merged_ranges = []
      
    for key, value in groupby(ranges, key_func):
        merged_ranges.extend(merge(list(value)))

    return merged_ranges


def get_hgvs_contextuals_url(hgvs):
    return f"{NCBI_VARIATION_SERVICES_BASE_URL}hgvs/{hgvs}/contextuals"


def get_spdi_all_equivalent_contextual_url(contextual_SPDI):
    return f'{NCBI_VARIATION_SERVICES_BASE_URL}spdi/{contextual_SPDI}/all_equivalent_contextual'


def get_spdi_canonical_representative_url(contextual_SPDI):
    return f'{NCBI_VARIATION_SERVICES_BASE_URL}spdi/{contextual_SPDI}/canonical_representative'


def build_spdi(seq_id, position, deleted_sequence, inserted_sequence):
    return f"{seq_id}:{position}:{deleted_sequence}:{inserted_sequence}"


def get_spdi_elements(response_object):
    return (response_object['seq_id'], response_object['position'], response_object['deleted_sequence'], response_object['inserted_sequence'])


def validate_subject(patient_id):
    if not patients_db.find_one({"patientID": patient_id}):
        abort(400, f"Patient ({patient_id}) not found.")


def get_variant(variant):
    variant = variant.strip()
    if variant.count(":") not in [1, 3]:
        abort(400, f'variant ({variant}) is not in the correct format(SPDI|HGVS)')

    variant = variant.lstrip()

    if variant.count(":") == 1: # HGVS expression
        SPDIs = hgvs_2_contextual_SPDIs(variant)
        if not SPDIs:
            abort(400, f'Cannot normalize variant: {variant}')
        elif not SPDIs["GRCh37"] and not SPDIs["GRCh38"]:
            abort(400, f'Cannot normalize variant: {variant}')
        else:
            normalized_variant = {"variant": variant, "GRCh37": SPDIs["GRCh37"], "GRCh38": SPDIs["GRCh38"]}

    elif variant.count(":") == 3: # SPDI expression
        SPDIs = SPDI_2_contextual_SPDIs(variant)
        if not SPDIs:
            abort(400, f'Cannot normalize variant: {variant}')
        elif not SPDIs["GRCh37"] and not SPDIs["GRCh38"]:
            abort(400, f'Cannot normalize variant: {variant}')
        else:
            normalized_variant = {"variant": variant, "GRCh37": SPDIs["GRCh37"], "GRCh38": SPDIs["GRCh38"]}
    else:
        abort(400, f'variant ({variant}) is not in the correct format(SPDI|HGVS)')

    return normalized_variant


def get_gene(gene):
    gene = gene.strip()
    gene_return = {'isSystem': False, 'gene': gene, 'system': None}
    if "|" in gene:
        if gene.count("|") == 1 and (re.match(SUPPORTED_GENE_SYSTEM_URLS, gene.rsplit('|')[0])):
            gene_system_url = gene.rsplit("|")[0]
            gene = gene.rsplit("|")[1]
            gene_return['isSystem'] = True
            gene_return['gene'] = gene
            gene_return['system'] = gene_system_url
        else:
            abort(400, f'gene ({gene}) is not in the correct format(codesystem|code)')

    return gene_return


def get_haplotype(haplotype):
    haplotype = haplotype.strip()
    haplotype_return = {'isSystem': False, 'haplotype': haplotype, 'system': None}
    if "|" in haplotype:
        if haplotype.count("|") == 1:
            haplotype_system_url = haplotype.rsplit("|")[0]
            haplotype = haplotype.rsplit("|")[1]
            haplotype_return['isSystem'] = True
            haplotype_return['haplotype'] = haplotype
            haplotype_return['system'] = haplotype_system_url
        else:
            abort(400, f'haplotype ({haplotype}) is not in the correct format(codesystem|code)')

    return haplotype_return


def get_treatment(treatment):
    treatment = treatment.strip()
    treatment_return = {'isSystem': False, 'treatment': (treatment if not is_int(treatment) else int(treatment)), 'system': None}
    if "|" in treatment:
        if treatment.count("|") == 1:
            treatment_system_url = treatment.rsplit("|")[0]
            treatment = int(treatment.rsplit("|")[1])
            treatment_return['isSystem'] = True
            treatment_return['treatment'] = treatment
            treatment_return['system'] = treatment_system_url
        else:
            abort(400, f'treatment ({treatment}) is not in the correct format(codesystem|code)')

    return treatment_return


def get_condition(condition):
    condition = condition.strip()
    condition_return = {'isSystem': False, 'condition': condition, 'system': None}
    if "|" in condition:
        if condition.count("|") == 1:
            condition_system_url = condition.rsplit("|")[0]
            condition = condition.rsplit("|")[1]
            condition_return['isSystem'] = True
            condition_return['condition'] = condition
            condition_return['system'] = condition_system_url
        else:
            abort(400, f'condition ({condition}) is not in the correct format(codesystem|code)')

    return condition_return


def get_range(range):
    range = range.strip()

    ref_seq = range.split(':')[0]
    chrom = get_build_and_chrom_by_ref_seq(ref_seq)

    if chrom is None:
        abort(400, f'RefSeq({ref_seq}) is not valid')

    chrom = chrom["chrom"]
    _range = list(map(int, range.split(':')[1].split('-')))

    if _range[0] > _range[1]:
        abort(400, f'Range start({_range[0]}) should be less than or equal to Range end({_range[1]})')

    chromosome = {'CHROM': chrom, 'RefSeq': ref_seq}
    _range = {'L': _range[0], 'H': _range[1]}

    return {'CHROMOSOME': chromosome, 'RANGE': _range}


def get_date(dateRange):
    dateRange = dateRange.strip()

    operator_code = dateRange[:2]

    operator = OPERATOR_CODE_TO_OPERATOR[operator_code]
    date_ = datetime.strptime(dateRange[2:], SUPPORTED_DATE_FORMAT)

    return {'OPERATOR': operator, 'DATE': date_}


def get_genomic_source_class(genomic_source_class):
    genomic_source_class = genomic_source_class.strip()

    if (genomic_source_class.lower() not in SUPPORTED_GENOMIC_SOURCE_CLASSES):
        abort(400, f"Genomic Source Class must be either 'germline' or 'somatic'")

    return genomic_source_class.lower()


def is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False


def get_genomics_build_presence(query):
    genomics_build_presence = {"GRCh37": False, "GRCh38": False}

    query["genomicBuild"] = "GRCh37"
    if variants_db.find_one(query):
        genomics_build_presence["GRCh37"] = True

    query["genomicBuild"] = "GRCh38"
    if variants_db.find_one(query):
        genomics_build_presence["GRCh38"] = True

    query.pop("genomicBuild")

    return genomics_build_presence


def get_genomics_build_presence_tests_db(query):
    genomics_build_presence = {"GRCh37": False, "GRCh38": False}

    query["genomicBuild"] = "GRCh37"
    if tests_db.find_one(query):
        genomics_build_presence["GRCh37"] = True

    query["genomicBuild"] = "GRCh38"
    if tests_db.find_one(query):
        genomics_build_presence["GRCh38"] = True

    query.pop("genomicBuild")

    return genomics_build_presence


def get_chromosome_to_ranges(ranges, genomics_build_presence):
    chromosome_to_ranges = []

    if ranges:
        for range_ in ranges:
            temp_dict = {}

            chromosome = range_['CHROMOSOME']
            _range = range_['RANGE']

            provided_genomic_build = get_build_and_chrom_by_ref_seq(chromosome['RefSeq'])["build"]
            other_genomic_build = get_other_build(provided_genomic_build)

            temp_dict["CHROM"] = chromosome['CHROM']
            temp_dict["RefSeq"] = chromosome['RefSeq']
            temp_dict["PGB"] = {"RefSeq": chromosome['RefSeq'], "BUILD": provided_genomic_build, "L": _range["L"], "H": _range["H"]}
            
            if (genomics_build_presence is None or genomics_build_presence[other_genomic_build]):
                rse_other_build = lift_over(chromosome['RefSeq'], _range["L"], _range["H"])

                if rse_other_build is None:
                    abort(422, f"Failed LiftOver ({chromosome['RefSeq']}:{_range['L']}-{_range['H']})")
                else:
                    other_ref_seq = get_ref_seq_by_chrom_and_build(other_genomic_build, chromosome['CHROM'])
                    temp_dict["OGB"] = {"RefSeq": other_ref_seq, "BUILD": other_genomic_build, "L": rse_other_build["start"], "H": rse_other_build["end"]}

            else:
                temp_dict["OGB"] = None

            chromosome_to_ranges.append(temp_dict)

    return chromosome_to_ranges


def get_dna_chg(svtype):
    dna_chg = SVTYPE_TO_DNA_CHANGE_TYPE.get(svtype.upper())
    return {"CODE": dna_chg[0], "DISPLAY": dna_chg[1]}


def create_fhir_variant_resource(record, subject):
    vid = f"dv-{str(record['_id'])}"

    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = vid
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://loinc.org",
                                   "code": "69548-6",
                                   "display": "Genetic variant assessment"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    resource["valueCodeableConcept"] = {"coding": [{"system": "http://loinc.org",
                                                    "code": "LA9633-4",
                                                    "display": "present"}]}

    resource["component"] = []

    if 'SVTYPE' in record:
        dna_chg = get_dna_chg(record['SVTYPE'])
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "48019-4",
                                                           "display": "DNA Change Type"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://sequenceontology.org",
                                                                           "code": f"{dna_chg['CODE']}",
                                                                           "display": f"{dna_chg['DISPLAY']}"}]}})


    # Genomic Source Class
    if 'genomicSourceClass' in record and record['genomicSourceClass'].upper() in ['GERMLINE', 'SOMATIC']:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "48002-0",
                                                           "display": "Genomic Source Class"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                           "code": f"{GENOMIC_SOURCE_CLASS_TO_CODE[record['genomicSourceClass'].upper()]}",
                                                                           "display": f"{record['genomicSourceClass']}"}]}})

    # Genomic reference sequence ID
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "48013-7",
                                                       "display": "Genomic reference sequence ID"}]},
                                  "valueCodeableConcept": {"coding": [{"system": "http://www.ncbi.nlm.nih.gov/nuccore",
                                                                       "code": f"{get_ref_seq_by_chrom_and_build(record['genomicBuild'], record['CHROM'])}"}]}})

    # Allelic State
    if "allelicState" in record and (('SVTYPE' not in record) or ('SVTYPE' in record and record['SVTYPE'] in ['DUP', 'DEL', 'INV', 'INS'] and record['genomicSourceClass'].lower() == 'germline')):
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "53034-5",
                                                           "display": "Allelic state"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                           "code": f"{ALLELIC_STATE_TO_CODE[record['allelicState'].upper()]}",
                                                                           "display": record["allelicState"]}]}})

    # Variation Code
    if "SPDI" in record:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81252-9",
                                                           "display": "Discrete genetic variant"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "https://api.ncbi.nlm.nih.gov/variation/v0/",
                                                                           "code": record["SPDI"],
                                                                           "display": record["SPDI"]}]}})
    
    # Allelic Frequency
    if ('SVTYPE' not in record):
        try:
            ad = float(record["ADS"][1]["AD"])
            dp = float(record["DP"])
            allelic_frequency = ad/dp
            resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                               "code": "81258-6",
                                                               "display": "Sample VAF"}]},
                                          "valueQuantity": {"value": allelic_frequency,
                                                            "unit": "relative frequency of a particular allele in the specimen",
                                                            "system": "http://unitsofmeasure.org",
                                                            "code": "1"}})
        except Exception as e:
            print(f"DEBUG: Error calculating Allelic Frequency - {e}")

    # Copy Number
    if (('SVTYPE' in record) and (record['SVTYPE'] in ['CNV', 'DUP', 'DEL']) and 'CN' in record):
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "82155-3",
                                                           "display": "Genomic Structural Variant copy Number"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://unitsofmeasure.org",
                                                                           "code": 1,
                                                                           "display": f"{record['CN']}"}]}})

    # Genomic Ref Allele ID
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                           "code": "69547-8",
                                           "display": "Genomic Ref allele [ID]"}]},
                                  "valueString": f"{record['REF']}"})


    # Genomic Alt Allele ID
    if ((('SVTYPE' not in record) or ('SVTYPE' in record and record['SVTYPE'] in ['INS'])) and ('ALT' in record) and (record['ALT'].isalpha())):
        alt = record['ALT']
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "69551-0",
                                                           "display": "Genomic Alt allele [ID]"}]},
                                      "valueString": f"{alt}"})

    # Genomic Coord System
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "92822-6",
                                                       "display": "Genomic coord system"}]},
                                  "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                       "code": "LA30100-4",
                                                                       "display": "0-based interval counting"}]}})

    # Variant Exact Start and End
    if 'SVTYPE' not in record:
        exact_start = record["POS"]
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81254-5",
                                                           "display": "Variant exact start-end"}]},
                                      "valueRange": {"low": {"value": exact_start}}})

    # Variant Outer/Inner Start and End
    if 'SVTYPE' in record:
        if 'CIPOS' in record and 'CIEND' in record:
            inner_start = record['POS'] + record['CIPOS'][1]
            inner_end = record['END'] - abs(record['CIEND'][0])
            outer_start = record['POS'] - abs(record['CIPOS'][0])
            outer_end = record['END'] + record['CIEND'][1]

            resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                               "code": "81301-4",
                                                               "display": "Variant outer start-end"}]},
                                          "valueRange": {"low": {"value": outer_start},
                                                         "high": {"value": outer_end}}})
        else:
            inner_start = record['POS']
            inner_end = record['END']

        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81302-2",
                                                           "display": "Variant inner start-end"}]},
                                      "valueRange": {"low": {"value": inner_start},
                                                     "high": {"value": inner_end}}})

    return resource


def create_dx_implication_profile(implication, subject, vids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(implication['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/diagnostic-implication"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "diagnostic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if len(vids) > 0:
        resource["derivedFrom"] = []
        for vid in vids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{vid}"})

    resource["component"] = []
    
    if implication['clinicalSignificance'] in CLIN_SIG_TO_CODE.keys():
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                            "code": "53037-8",
                                                            "display": "Genetic variation clinical significance"}]},
                                        "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                            "code": f"{CLIN_SIG_TO_CODE[implication['clinicalSignificance']]}",
                                                                            "display": f"{implication['clinicalSignificance']}"}]}})
    else:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                            "code": "53037-8",
                                                            "display": "Genetic variation clinical significance"}]},
                                        "valueCodeableConcept": {"text": f"{implication['clinicalSignificance']}"}})

    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81259-4",
                                                           "display": "predicted phenotype"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "https://www.ncbi.nlm.nih.gov/medgen",
                                                                           "code": f"{implication['MedGenCode']}",
                                                                           "display": f"{implication['MedGenDesc']}"}]}})

    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "93044-6",
                                                           "display": "Level of evidence"}]},
                                      "valueCodeableConcept": {"text": f"{implication['evidenceLevel']}"}})
    
    return resource


def create_tx_implication_profile_civic(implication, subject, vids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(implication['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/therapeutic-implication"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "therapeutic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if len(vids) > 0:
        resource["derivedFrom"] = []
        for vid in vids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{vid}"})


    resource["component"] = []

    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "81259-4",
                                                       "display": "phenotypic treatment context"}]},
                                           "valueCodeableConcept": {"coding": [{"system": "https://disease-ontology.org",
                                                                                "code": f"{implication['diseaseOntologyCode']}",
                                                                                "display": f"{implication['diseaseOntologyDesc']}"}]}})

    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "93044-6",
                                                       "display": "Level of evidence"}]},
                                           "valueCodeableConcept": {"text": f"{implication['evidenceLevel']}"}})
    
    for med in implication['medicationAssessed']:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "51963-7",
                                                           "display": "medication-assessed"}]},
                    "valueCodeableConcept": {"coding": [{"system": "http://www.nlm.nih.gov/research/umls/rxnorm",
                                                         "code": f"{med['rxnormCode']}",
                                                         "display": f"{med['rxnormDesc']}"}]}})
        
    resource["component"].append({"code": {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                                       "code": "predicted-therapeutic-implication",
                                                       "display": "predicted-therapeutic-implication"}]},
                "valueCodeableConcept": {"text": implication['predictedImplication']}})
    
    return resource


def create_tx_implication_profile_pharmgkb(implication, subject, gids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(implication['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/therapeutic-implication"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "therapeutic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if len(gids) > 0:
        resource["derivedFrom"] = []
        for gid in gids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{gid}"})

    resource["component"] = []

    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "93044-6",
                                                       "display": "Level of evidence"}]},
                                           "valueCodeableConcept": {"text": f"{implication['evidenceLevel']}"}})
    
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                        "code": "51963-7",
                                                        "display": "medication-assessed"}]},
                "valueCodeableConcept": {"coding": [{"system": "http://www.nlm.nih.gov/research/umls/rxnorm",
                                                        "code": f"{implication['rxnormCode']}",
                                                        "display": f"{implication['rxnormDesc']}"}]}})
        
    resource["component"].append({"code": {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                                       "code": "predicted-therapeutic-implication",
                                                       "display": "predicted-therapeutic-implication"}]},
                "valueCodeableConcept": {"text": implication['predictedImplication']}})
    
    return resource


def create_haplotype_profile(haplotype, subject, gid):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(haplotype['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/haplotype"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "diagnostic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if gid:
        resource["derivedFrom"] = [{"reference": f"Observation/{gid}"}]

    return resource


def create_genotype_profile(genotype, subject, gids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(genotype['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genotype"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://loinc.org",
                                                       "code": "84413-4",
                                                       "display": "Genotype display name"}]}

    resource["valueCodeableConcept"] = {"coding": [{"system": f"{genotype['genotypeCodeSystem']}",
                                         "code": f"{genotype['genotypeCode']}",
                                         "display": f"{genotype['genotypeDesc']}"}]}

    if len(gids) > 0:
        resource["derivedFrom"] = []
        for gid in gids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{gid}"})

    resource["component"] = []
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "48018-6",
                                                           "display": "Gene studied"}]},

                                  "valueCodeableConcept": {"coding": [{"system": "http://www.genenames.org/geneId",
                                                                        "code": f"{genotype['geneCode']}",
                                                                        "display": f"{genotype['geneDesc']}"}]}})
    
    return resource


def add_variation_id(resource, variation_id, system):
    spdi_index = next((i for i, item in enumerate(resource["component"]) if item["code"]["coding"][0]["code"] == "81252-9"), None)

    if spdi_index is not None:
        resource["component"][spdi_index]["valueCodeableConcept"]["coding"].append(
            {"system": system,
             "code": f'{variation_id}'})
    else:
        resource["component"].append(
            {"code": {"coding": [{"system": "http://loinc.org",
                                  "code": "81252-9",
                                  "display": "Discrete genetic variant"}]},
             "valueCodeableConcept": {"coding": [{"system": system,
                                                  "code": f'{variation_id}'}]}})


def lift_over(ref_seq, start, end):
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        return None

    tempDict = get_build_and_chrom_by_ref_seq(ref_seq)
    if tempDict:
        chrom = tempDict['chrom']
        build = tempDict['build']
    else:
        return None

    # convert build37 to build38
    if build == 'GRCh37':
        new_ref_seq = get_ref_seq_by_chrom_and_build('GRCh38', chrom)
        lo_b37_to_b38 = pyliftover.LiftOver('hg19', 'hg38')
        new_coor = lo_b37_to_b38.convert_coordinate(chrom, start, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_start = new_coor[0][1]            
        new_coor = lo_b37_to_b38.convert_coordinate(chrom, end, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_end = new_coor[0][1]      
        if new_start > new_end:
            return None
        else:
            return {"refSeq": new_ref_seq, "start": new_start, "end": new_end}

    #convert build38 to build37
    elif build == 'GRCh38':
        new_ref_seq = get_ref_seq_by_chrom_and_build('GRCh37', chrom)
        lo_b38_to_b37 = pyliftover.LiftOver('hg38', 'hg19')
        
        new_coor = lo_b38_to_b37.convert_coordinate(chrom, start, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_start = new_coor[0][1]    
        new_coor = lo_b38_to_b37.convert_coordinate(chrom, end, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_end = new_coor[0][1]
        if new_start > new_end:
            return None
        else:
            return {"refSeq": new_ref_seq, "start":new_start, "end": new_end}
    else:
        return None


def get_build_and_chrom_by_ref_seq(ref_seq):
    chrom_data = chromosomes_data.find_one({'$or': [{'build37RefSeq': {'$eq': ref_seq}}, {'build38RefSeq': {'$eq': ref_seq}}]})

    if chrom_data:
        if chrom_data['build37RefSeq'] == ref_seq:
            return {"chrom": chrom_data['chr'], "build": 'GRCh37'}
        else:
            return {"chrom": chrom_data['chr'], "build": 'GRCh38'}
    
    return None


def get_ref_seq_by_chrom_and_build(build, chrom):
    chrom_data = chromosomes_data.find_one({'chr': {'$eq': chrom}})

    if chrom_data:
        if build == 'GRCh37':
            return chrom_data['build37RefSeq']
        if build == 'GRCh38':
            return chrom_data['build38RefSeq']

    return None


def get_other_build(build):
    return "GRCh37" if build == 'GRCh38' else "GRCh38"


def get_intersected_regions(bed_id, build, chrom, start, end, intersected_regions):
    try:
        result = beds.aggregate([
            {"$match": {"$and": [{"BedID": {"$eq": bed_id}},
                                 {"BED.Chromosome": {"$eq": chrom}},
                                 {'$or': [
                                    {"$and": [
                                        {"BED.Start": {"$gte": start}},
                                        {"BED.Start": {"$lte": end}}
                                    ]},
                                    {"$and": [
                                        {"BED.Start": {"$lte": start}},
                                        {"BED.End": {"$gte": start}}
                                    ]}]}]}},
            {"$unwind": "$BED"},
            {"$match": {"$and": [{"BED.Chromosome": {"$eq": chrom}},
                                 {'$or': [
                                    {"$and": [
                                        {"BED.Start": {"$gte": start}},
                                        {"BED.Start": {"$lte": end}}
                                    ]},
                                    {"$and": [
                                        {"BED.Start": {"$lte": start}},
                                        {"BED.End": {"$gte": start}}
                                    ]}]}]}},
            {"$group": {"_id": "$BedID", "BED": {"$push": "$$ROOT.BED"}}}])

        result = list(result)
    except Exception as e:
        print(f"DEBUG: Error under get_intersected_regions(bed_id={bed_id}, build={build}, chrom={chrom}, start={start}, end={end}, intersected_regions={intersected_regions})")
        result = result

    if result:
        result = result[0]["BED"]

        for csePair in result:
            ref_seq = get_ref_seq_by_chrom_and_build(build, csePair["Chromosome"])
            intersected_regions.append(f'{ref_seq}:{max(start, csePair["Start"])}-{min(end, csePair["End"])}')


def hgvs_2_contextual_SPDIs(hgvs):

    # convert hgvs to contextualSPDI
    url = get_hgvs_contextuals_url(hgvs)
    headers = {'Accept': 'application/json'}

    r = requests.get(url, headers=headers)
    if r.status_code != 200: 
        return False

    response = r.json()
    raw_data = response['data']
    raw_SPDI = raw_data['spdis'][0]

    seq_id, position, deleted_sequence, inserted_sequence = get_spdi_elements(raw_SPDI)

    contextual_SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)
    
    # convert contextualSPDI to build37 and build38 contextual SPDIs
    url = get_spdi_all_equivalent_contextual_url(contextual_SPDI)
    headers = {'Accept': 'application/json'}

    r = requests.get(url, headers=headers)
    if r.status_code != 200: 
        return False

    response = r.json()
    raw_SPDI_List = response['data']['spdis']

    b37SPDI = None
    b38SPDI = None
    for item in raw_SPDI_List:
        if item['seq_id'].startswith("NC_"):
            temp = get_build_and_chrom_by_ref_seq(item['seq_id'])
            if temp:
                seq_id, position, deleted_sequence, inserted_sequence = get_spdi_elements(item)

                if temp['build']=='GRCh37':
                    b37SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)
                elif temp['build'] == 'GRCh38':
                    b38SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)
            else:
                return False

    return {"GRCh37": b37SPDI, "GRCh38": b38SPDI}


def hgvs_2_canonical_SPDI(hgvs):

    # convert hgvs to contextualSPDI
    url = get_hgvs_contextuals_url(hgvs)
    headers = {'Accept': 'application/json'}

    r = requests.get(url, headers=headers)
    if r.status_code != 200: 
        return False

    response = r.json()
    raw_data = response['data']
    raw_SPDI = raw_data['spdis'][0]

    seq_id, position, deleted_sequence, inserted_sequence = get_spdi_elements(raw_SPDI)

    contextual_SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)
    
    # convert contextualSPDI to canonical SPDI
    url = get_spdi_canonical_representative_url(contextual_SPDI)
    headers = {'Accept': 'application/json'}

    r = requests.get(url, headers=headers)
    if r.status_code != 200: 
        return False

    response = r.json()
    raw_SPDI = response['data']

    seq_id, position, deleted_sequence, inserted_sequence = get_spdi_elements(raw_SPDI)

    canonical_SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)

    return {"canonicalSPDI": canonical_SPDI}


def SPDI_2_contextual_SPDIs(spdi):
    url = get_spdi_all_equivalent_contextual_url(spdi)
    headers = {'Accept': 'application/json'}

    r = requests.get(url, headers=headers)
    if r.status_code != 200:
        return False

    response = r.json()
    raw_SPDI_List = response['data']['spdis']

    b37SPDI = None
    b38SPDI = None
    for item in raw_SPDI_List:
        if item['seq_id'].startswith("NC_"):
            temp = get_build_and_chrom_by_ref_seq(item['seq_id'])
            if temp:
                seq_id, position, deleted_sequence, inserted_sequence = get_spdi_elements(item)

                if temp['build'] == 'GRCh37':
                    b37SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)
                elif temp['build'] == 'GRCh38':
                    b38SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)
            else:
                return False

    return {"GRCh37": b37SPDI, "GRCh38": b38SPDI}


def SPDI_2_canonical_SPDI(spdi):
    url = get_spdi_canonical_representative_url(spdi)
    headers = {'Accept': 'application/json'}

    r = requests.get(url, headers=headers)
    if r.status_code != 200: 
        return False

    response = r.json()
    raw_SPDI = response['data']

    seq_id, position, deleted_sequence, inserted_sequence = get_spdi_elements(raw_SPDI)

    canonical_SPDI = build_spdi(seq_id, position, deleted_sequence, inserted_sequence)

    return {"canonicalSPDI": canonical_SPDI}


def query_clinvar_by_variants(normalized_variant_list, code_list, query, population=False):
    variant_list = []
    for item in normalized_variant_list:
        if item["GRCh37"]:
            variant_list.append(item["GRCh37"])
        if item["GRCh38"]:
            variant_list.append(item["GRCh38"])

    pipeline_part = [{'$match': {'$expr': {'$and': [{'$or': [{'$eq': ['$b37SPDI', '$$mySPDI']},
                                                             {'$eq': ['$b38SPDI', '$$mySPDI']}]}]}}},
                     {'$addFields': {}}]
    if code_list != []:
        pipeline_part.append({'$match': {'$or': []}})
        for condition in code_list:
            if condition['isSystem']:
                pipeline_part[-1]['$match']['$or'].append({'MedGenCode': {'$eq': condition['condition']}})
            else:
                pipeline_part[-1]['$match']['$or'].append({'$or': [
                    {'MedGenCode': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                    {'MedGenDesc': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                ]})

    query['SPDI'] = {'$in': variant_list}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'Clinvar', 'let': {'mySPDI': '$SPDI'}, 'pipeline': pipeline_part,
                                 'as': 'clinvarMatches'}},
                    {'$addFields': {}},
                    {'$match': {'clinvarMatches': {'$exists': True, '$not': {'$size': 0}}}}]
             
    if population:
        query_string.append({'$group': {'_id': '$patientID'}})

    try:
        results = variants_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_clinvar_by_variants(normalized_variant_list={normalized_variant_list}, code_list={code_list}, query={query}, population={population})")
        results = []

    query_results = []
    
    for item in results:
        item["UUID"] = str(uuid4())
        query_results.append(item)
  
    return query_results


def query_clinvar_by_condition(code_list, query):
    condition_query = {'$or': []}

    for condition in code_list:
        if condition['isSystem']:
            condition_query['$or'].append({'MedGenCode': {'$eq': condition['condition']}})
        else:
            condition_query['$or'].append({'$or': [
                {'MedGenCode': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                {'MedGenDesc': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
            ]})

    query_string = [
        {'$match': condition_query},
        {'$lookup': {'from': 'Variants', 'let': {'b37SPDI': '$b37SPDI', 'b38SPDI': '$b38SPDI'}, 
                     'pipeline': [{'$match': query},
                                  {'$addFields': {}},
                                  {'$match': {'$expr': {'$and': [{'$or': [
                                                                    {'$eq': ['$SPDI', '$$b37SPDI']},
                                                                    {'$eq': ['$SPDI', '$$b38SPDI']}
                                                                ]}]}}},
                                  {'$addFields': {}}],
                      'as': 'patientMatches'}},
        {'$addFields': {}}, 
        {'$match': {'patientMatches': {'$exists': True, '$not': {'$size': 0}}}}
    ]
 
    try:
        results = clinvar_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_clinvar_by_condition(code_list={code_list}, query={query})")
        results = []

    query_results = []

    for item in results:
        for match in item["patientMatches"]:
            match["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_CIVIC_by_variants(normalized_variant_list, code_list, treatment_list, query, population=False):
    variant_list = []
    for item in normalized_variant_list:
        if item["GRCh37"]:
            variant_list.append(item["GRCh37"])
        if item["GRCh38"]:
            variant_list.append(item["GRCh38"])

    pipeline_part = [{'$match': {'$expr': {'$and': [{'$or': [{'$eq': ['$b37SPDI', '$$mySPDI']},
                                                             {'$eq': ['$b38SPDI', '$$mySPDI']}]}]}}},
                     {'$addFields': {}}]

    if code_list != []:
        pipeline_part.append({'$match': {'$or': []}})
        for condition in code_list:
            if condition['isSystem']:
                pipeline_part[-1]['$match']['$or'].append({'diseaseOntologyCode': {'$eq': condition['condition']}})
            else:
                pipeline_part[-1]['$match']['$or'].append({'$or': [
                    {'diseaseOntologyCode': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                    {'diseaseOntologyDesc': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                ]})

    if treatment_list != []:
        pipeline_part.append({'$match': {'$or': []}})
        for treatment in treatment_list:
            if treatment['isSystem']:
                pipeline_part[-1]['$match']['$or'].append({'medicationAssessed.rxnormCode': {'$eq': f"{treatment['treatment']}"}})
            else:
                pipeline_part[-1]['$match']['$or'].append({'$or': [
                    {'medicationAssessed.rxnormCode': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                    {'medicationAssessed.rxnormDesc': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                ]})

    query['SPDI'] = {'$in': variant_list}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'CIViC', 'let': {'mySPDI': '$SPDI'}, 'pipeline': pipeline_part,
                                 'as': 'civicMatches'}},
                    {'$addFields': {}}, 
                    {'$match': {'civicMatches': {'$exists': True, '$not': {'$size': 0}}}}]
             
    if population:
        query_string.append({'$group': {'_id': '$patientID'}})

    try:
        results = variants_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_CIVIC_by_variants(normalized_variant_list={normalized_variant_list}, code_list={code_list}, treatment_list={treatment_list}, query={query}, population={population})")
        results = []

    query_results = []
    
    for item in results:
        item["UUID"] = str(uuid4())
        query_results.append(item)
  
    return query_results


def query_CIVIC_by_condition(code_list, treatment_list, query):
    match_dict = {"$and": []}
    if code_list != []:
        match_dict['$and'].append({'$or': []})
        for condition in code_list:
            if condition['isSystem']:
                match_dict['$and'][0]['$or'].append({'diseaseOntologyCode': {'$eq': condition['condition']}})
            else:
                match_dict['$and'][0]['$or'].append({'$or': [
                    {'diseaseOntologyCode': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                    {'diseaseOntologyDesc': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                ]})

    if treatment_list != []:
        match_dict['$and'].append({'$or': []})
        for treatment in treatment_list:
            if treatment['isSystem']:
                match_dict['$and'][-1]['$or'].append({'medicationAssessed.rxnormCode': {'$eq': f"{treatment['treatment']}"}})
            else:
                match_dict['$and'][-1]['$or'].append({'$or': [
                    {'medicationAssessed.rxnormCode': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                    {'medicationAssessed.rxnormDesc': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                ]})

    query_string = [{'$match': match_dict},
                    {'$lookup': {'from': 'Variants', 'let': {'b37SPDI': '$b37SPDI', 'b38SPDI': '$b38SPDI'}, 
                                 'pipeline': [{'$match': query},
                                              {'$match': {'$expr': {'$and': [{'$or': [
                                                                               {'$eq': ['$SPDI', '$$b37SPDI']},
                                                                               {'$eq': ['$SPDI', '$$b38SPDI']}
                                                                           ]}]}}},
                                              {'$addFields': {}}],
                                 'as': 'patientMatches'}},
                    {'$addFields': {}},
                    {'$match': {'patientMatches': {'$exists': True, '$not': {'$size': 0}}}}]

    try:
        results = civic_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_CIVIC_by_condition(code_list={code_list}, treatment_list={treatment_list}, query={query})")
        results = []

    query_results = []

    for item in results:
        for match in item["patientMatches"]:
            match["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_PharmGKB_by_haplotypes(normalizedHaplotypeList, treatmentCodeList, query, population=False):
    pipeline_part = [{'$match': {'$expr': {'$and': [{'$eq': ['$genotypeCode', '$$mygenotypeCode']}]}}},
                    {'$addFields': {}}]

    if treatmentCodeList != []:
        pipeline_part.append({'$match': {'$or': []}})
        for treatment in treatmentCodeList:
            if treatment['isSystem']:
                pipeline_part[-1]['$match']['$or'].append({'rxnormCode': {'$eq': treatment['treatment']}})
            else:
                pipeline_part[-1]['$match']['$or'].append({'$or': [
                        {'rxnormCode': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                        {'rxnormDesc': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                    ]})
    
    query['$or'] = []

    for haplotype in normalizedHaplotypeList:

        if haplotype['isSystem']:
            query['$or'].append({'genotypeCode': {"$eq": haplotype['haplotype']}})
        else:
            query['$or'].append({'$or': [
                    {'genotypeCode': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}},
                    {'genotypeDesc': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}}
                ]})

    query_string = [{'$match': query},
                   {'$lookup': {'from': 'PharmGKB', 'let': {'mygenotypeCode': '$genotypeCode'}, 'pipeline': pipeline_part,
                                'as': 'pahrmGKBMatches'}},
                   {'$addFields': {}},
                   {'$match': {'pahrmGKBMatches': {'$exists': True, '$not': {'$size': 0}}}}]

    if population:
        query_string.append({'$group': {'_id': '$patientID'}})
             
    try:
        results = genotypes_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_PharmGKB_by_haplotypes(normalizedHaplotypeList={normalizedHaplotypeList}, treatment_list={treatment_list}, query={query}, population={population})")
        results = []

    query_results = []
    
    for item in results:
        item["UUID"] = str(uuid4())
        query_results.append(item)
  
    return query_results


def query_PharmGKB_by_treatments(code_list, treatment_list, query):
    match_dict = {"$and": []}
    if code_list != []:
        match_dict['$and'].append({'$or': []})
        for condition in code_list:
            if condition['isSystem']:
                match_dict['$and'][0]['$or'].append({'diseaseOntologyCode': {'$eq': condition['condition']}})
            else:
                match_dict['$and'][0]['$or'].append({'$or': [
                        {'diseaseOntologyCode': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                        {'diseaseOntologyDesc': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                    ]})

    if treatment_list != []:
        match_dict['$and'].append({'$or': []})
        for treatment in treatment_list:
            if treatment['isSystem']:
                match_dict['$and'][-1]['$or'].append({'rxnormCode': {'$eq': treatment['treatment']}})
            else:
                match_dict['$and'][-1]['$or'].append({'$or': [
                        {'rxnormCode': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                        {'rxnormDesc': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                    ]})
    
    if 'genomicSourceClass' in query:
        query.pop('genomicSourceClass')
 
    query_string = [{'$match': match_dict},
                    {'$lookup': {'from': 'Genotypes', 'let': {'genotypeCode': '$genotypeCode'}, 
                                 'pipeline': [{'$match': query},
                                              {'$match': {'$expr': {'$and': [{'$eq': ['$genotypeCode', '$$genotypeCode']}]}}},
                                              {'$addFields': {}}],
                                 'as': 'patientMatches'}},
                    {'$addFields': {}},
                    {'$match': {'patientMatches': {'$exists': True, '$not': {'$size': 0}}}}]

    try:
        results = pharmGKB_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_PharmGKB_by_treatments(code_list={code_list}, treatment_list={treatment_list}, query={query})")
        results = []

    query_results = []

    for item in results:
        for match in item["patientMatches"]:
            match["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_genes(gene):
    query = {}

    pipeline_part = [{'$match': {'$expr': {'$eq': ['$ncbiGeneSymbol', '$$geneSymbol']}}},
                     {'$addFields': {}}]

    query['$or'] = [{"HGNCgeneId": {"$eq": gene}},
                    {"HGNCgeneSymbol": {"$eq": gene}}]

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'Transcripts', 'let': {'geneSymbol': '$NCBIgeneSymbol'}, 'pipeline': pipeline_part,
                                 'as': 'transcriptMatches'}},
                    {'$addFields': {}}]

    try:
        results = genes_data.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under get_feature_coordinates(gene={gene})")
        result = []
  
    return results


def query_transcript(transcript):
    query = {}

    gene_pipeline_part = [{'$match': {'$expr': {'$eq': ['$NCBIgeneSymbol', '$$geneSymbol']}}},
                          {'$addFields': {}}]

    exon_pipeline_part = [{'$match': {'$expr': {'$eq': ['$transcript', '$$refSeq']}}},
                          {'$addFields': {}}]

    query['transcriptRefSeq'] = {"$regex": ".*"+str(transcript).replace('*', r'\*')+".*"}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'Genes', 'let': {'geneSymbol': '$ncbiGeneSymbol'}, 'pipeline': gene_pipeline_part,
                                 'as': 'geneMatches'}},
                    {'$addFields': {}},
                    {'$lookup': {'from': 'Exons', 'let': {'refSeq': '$transcriptRefSeq'}, 'pipeline': exon_pipeline_part,
                                 'as': 'exonMatches'}},
                    {'$addFields': {}}]

    try:
        results = transcripts_data.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under get_feature_coordinates(transcript={transcript})")
        result = []

    if len(results) > 1:
        abort(400, "Unable to provide information on this transcript at this time")
  
    return results
