import json
from collections import OrderedDict

import requests
from flask import abort, jsonify

from app import common, input_normalization
from utilities import hla
from utilities.spdi_refseq import get_ref_seq_subseq


def fetch_concept_map(mapID):
    url = f"http://hapi.fhir.org/baseR4/ConceptMap/{mapID}"
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)


def SPDI_all_equivalent_contextual(variant):
    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/spdi/' + variant + '/all_equivalent_contextual'
    headers = {'Accept': 'application/json'}
    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        abort(500, "NCBI Variation Services server is down")
    return response


def SPDI_hgvs(variant):
    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/spdi/' + variant + '/hgvs'
    headers = {'Accept': 'application/json'}
    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        abort(500, "NCBI Variation Services issue")
    return response


def HGVS_contextuals(variant):
    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/hgvs/' + variant + '/contextuals'
    headers = {'Accept': 'application/json'}
    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        abort(500, "NCBI Variation Services issue")
    return response


def get_feature_coordinates(
        chromosome=None, gene=None, transcript=None, protein=None):

    if not chromosome and not gene and not transcript and not protein:
        abort(400, "You must provide one parameter.")

    if chromosome:
        if gene or transcript or protein:
            abort(400, "Only one argument can be suuplied")

        chromosome = chromosome.strip()

        try:
            result = common.chromosomes_data.aggregate([{"$match": {"chr": {"$eq": f"chr{chromosome[3:].upper()}"}}}])
            result = list(result)
        except Exception as e:
            print(f"DEBUG: Error({e}) under get_feature_coordinates(chromosome={chromosome})")
            result = []

        if not result:
            return jsonify([])

        result = result[0]

        output = OrderedDict()
        output["chromosome"] = result["chr"]
        output["chromosomeLink"] = f'https://medlineplus.gov/genetics/chromosome/{chromosome[3:].lower() if chromosome[3:].lower() != "m" else "mitochondrial-dna"}/'
        output["build37RefSeq"] = result["build37RefSeq"]
        output["build38RefSeq"] = result["build38RefSeq"]

        return jsonify(output)

    if gene:
        if transcript or protein:
            abort(400, "Only one argument can be suuplied")

        gene = gene.strip().upper()

        if common.is_int(gene):
            gene = f"HGNC:{gene}"

        result = common.query_genes(gene)

        output = []

        for res in result:
            ord_dict = OrderedDict()

            ord_dict['geneId'] = res['HGNCgeneId']
            ord_dict['geneSymbol'] = res['HGNCgeneSymbol']
            ord_dict['geneLink'] = f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{res['HGNCgeneId']}"
            ord_dict['build37Coordinates'] = f"{res['build37RefSeq']}:{res['build37Start']}-{res['build37End']}"
            ord_dict['build38Coordinates'] = f"{res['build38RefSeq']}:{res['build38Start']}-{res['build38End']}"
            ord_dict['transcripts'] = []
            ord_dict['MANE'] = []
            for transcript in res['transcriptMatches']:
                ord_dict['transcripts'].append(transcript['transcriptRefSeq'])
                if transcript['MANE'] == 1:
                    ord_dict['MANE'].append(transcript['transcriptRefSeq'])

            if not ord_dict['transcripts']:
                ord_dict.pop('transcripts')
            if not ord_dict['MANE']:
                ord_dict.pop('MANE')

            output.append(ord_dict)

        return (jsonify(output))

    if transcript:
        if protein:
            abort(400, "Only one argument can be suuplied")

        transcript = transcript.strip().upper()

        if '.' in transcript:
            transcript = transcript.split('.')[0]

        result = common.query_transcript(transcript)

        output = []

        for res in result:
            ord_dict = OrderedDict()

            ord_dict['transcript'] = res['transcriptRefSeq']
            ord_dict['transcriptLink'] = f"https://www.ncbi.nlm.nih.gov/nuccore/{res['transcriptRefSeq']}"

            if res['geneMatches']:
                ord_dict['geneSymbol'] = res['geneMatches'][0]['HGNCgeneSymbol']

            ord_dict['build37Coordinates'] = f"{res['build37RefSeq']}:{res['build37Start']}-{res['build37End']}"
            ord_dict['build38Coordinates'] = f"{res['build38RefSeq']}:{res['build38Start']}-{res['build38End']}"
            ord_dict['MANE'] = True if res['MANE'] == 1 else False

            if res['cdsMatches']:
                ord_dict['build38CDSStart'] = res['cdsMatches'][0]['build38CDSStart']
                ord_dict['build38CDSEnd'] = res['cdsMatches'][0]['build38CDSEnd']
                ord_dict['build38Strand'] = res['cdsMatches'][0]['build38Strand']
                ord_dict['build37CDSStart'] = res['cdsMatches'][0]['build37CDSStart']
                ord_dict['build37CDSEnd'] = res['cdsMatches'][0]['build37CDSEnd']
                ord_dict['build37Strand'] = res['cdsMatches'][0]['build37Strand']

            ord_dict['exons'] = []
            for exon in res['exonMatches']:
                exon_dict = OrderedDict()

                exon_dict['exonNumber'] = exon['exonNumber']
                exon_dict['build37Coordinates'] = f"{exon['build37RefSeq']}:{exon['build37Start']}-{exon['build37End']}"
                exon_dict['build38Coordinates'] = f"{exon['build38RefSeq']}:{exon['build38Start']}-{exon['build38End']}"

                ord_dict['exons'].append(exon_dict)

            if not ord_dict['exons']:
                ord_dict.pop('exons')
            else:
                ord_dict['exons'] = sorted(ord_dict['exons'], key=lambda d: d['exonNumber'])

            output.append(ord_dict)

        return (jsonify(output))

    if protein:

        protein = protein.strip().upper()

        if '.' in protein:
            protein = protein.split('.')[0]

        try:
            result = common.proteins_data.aggregate(
                [{"$match": {"proteinRefSeq": {'$regex': ".*"+str(protein).replace('*', r'\*')+".*"}}}])
            result = list(result)
        except Exception as e:
            print(f"DEBUG: Error({e}) under get_feature_coordinates(protein={protein})")
            result = []

        if len(result) > 1:
            abort(400, "Unable to provide information on this transcript at this time")

        if not result:
            return jsonify([])

        result = result[0]

        output = OrderedDict()
        output["protein"] = result["proteinRefSeq"]
        output["proteinLink"] = f'https://www.ncbi.nlm.nih.gov/protein/{output["protein"]}/'
        output["proteinName"] = result["proteinName"]
        output["transcript"] = result["transcript"]

        return (jsonify(output))


def find_the_gene(range=None):
    if not range:
        abort(400, "Range is required.")

    given_range = common.get_range(range)

    result = common.query_genes_range(given_range)

    output = []

    for res in result:
        ord_dict = OrderedDict()

        ord_dict['geneId'] = res['HGNCgeneId']
        ord_dict['geneSymbol'] = res['HGNCgeneSymbol']
        ord_dict['geneLink'] = f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{res['HGNCgeneId']}"
        ord_dict['build37Coordinates'] = f"{res['build37RefSeq']}:{res['build37Start']}-{res['build37End']}"
        ord_dict['build38Coordinates'] = f"{res['build38RefSeq']}:{res['build38Start']}-{res['build38End']}"
        ord_dict['transcripts'] = []
        ord_dict['MANE'] = []
        for transcript in res['transcriptMatches']:
            ord_dict['transcripts'].append(transcript['transcriptRefSeq'])
            if transcript['MANE'] == 1:
                ord_dict['MANE'].append(transcript['transcriptRefSeq'])

        if not ord_dict['transcripts']:
            ord_dict.pop('transcripts')
        if not ord_dict['MANE']:
            ord_dict.pop('MANE')

        output.append(ord_dict)

    return (jsonify(output))


def normalize_variant(variant):
    normalizedVariant = {
        "b37SPDI": "",
        "b38SPDI": "",
        "b37HGVS": "",
        "b38HGVS": ""
    }
    if variant.count(':') == 3:  # SPDI
        response = SPDI_all_equivalent_contextual(variant)
        for spdi in json.loads(response.text)["data"]["spdis"]:
            if spdi["seq_id"] in ['NC_000001.10', 'NC_000002.11', 'NC_000003.11', 'NC_000004.11', 'NC_000005.9', 'NC_000006.11', 'NC_000007.13', 'NC_000008.10', 'NC_000009.11', 'NC_000010.10', 'NC_000011.9', 'NC_000012.11', 'NC_000013.10', 'NC_000014.8', 'NC_000015.9', 'NC_000016.9', 'NC_000017.10', 'NC_000018.9', 'NC_000019.9', 'NC_000020.10', 'NC_000021.8', 'NC_000022.10', 'NC_000023.10', 'NC_000024.9']:
                normalizedVariant["b37SPDI"] = spdi["seq_id"] + ":" + str(spdi["position"]) + ":" + spdi["deleted_sequence"] + ":" + spdi["inserted_sequence"]
                normalizedVariant["b37HGVS"] = json.loads(SPDI_hgvs(normalizedVariant["b37SPDI"]).text)["data"]["hgvs"]
            if spdi["seq_id"] in ['NC_000001.11', 'NC_000002.12', 'NC_000003.12', 'NC_000004.12', 'NC_000005.10', 'NC_000006.12', 'NC_000007.14', 'NC_000008.11', 'NC_000009.12', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12', 'NC_000013.11', 'NC_000014.9', 'NC_000015.10', 'NC_000016.10', 'NC_000017.11', 'NC_000018.10', 'NC_000019.10', 'NC_000020.11', 'NC_000021.9', 'NC_000022.11', 'NC_000023.11', 'NC_000024.10']:
                normalizedVariant["b38SPDI"] = spdi["seq_id"] + ":" + str(spdi["position"]) + ":" + spdi["deleted_sequence"] + ":" + spdi["inserted_sequence"]
                normalizedVariant["b38HGVS"] = json.loads(SPDI_hgvs(normalizedVariant["b38SPDI"]).text)["data"]["hgvs"]
    elif variant.count(":") == 1:  # HGVS
        spdi = json.loads(HGVS_contextuals(variant).text)["data"]["spdis"][0]
        spdi = spdi["seq_id"] + ":" + str(spdi["position"]) + ":" + spdi["deleted_sequence"] + ":" + spdi["inserted_sequence"]
        response = SPDI_all_equivalent_contextual(spdi)
        for spdi in json.loads(response.text)["data"]["spdis"]:
            if spdi["seq_id"] in ['NC_000001.10', 'NC_000002.11', 'NC_000003.11', 'NC_000004.11', 'NC_000005.9', 'NC_000006.11', 'NC_000007.13', 'NC_000008.10', 'NC_000009.11', 'NC_000010.10', 'NC_000011.9', 'NC_000012.11', 'NC_000013.10', 'NC_000014.8', 'NC_000015.9', 'NC_000016.9', 'NC_000017.10', 'NC_000018.9', 'NC_000019.9', 'NC_000020.10', 'NC_000021.8', 'NC_000022.10', 'NC_000023.10', 'NC_000024.9']:
                normalizedVariant["b37SPDI"] = spdi["seq_id"] + ":" + str(spdi["position"]) + ":" + spdi["deleted_sequence"] + ":" + spdi["inserted_sequence"]
                normalizedVariant["b37HGVS"] = json.loads(SPDI_hgvs(normalizedVariant["b37SPDI"]).text)["data"]["hgvs"]
            if spdi["seq_id"] in ['NC_000001.11', 'NC_000002.12', 'NC_000003.12', 'NC_000004.12', 'NC_000005.10', 'NC_000006.12', 'NC_000007.14', 'NC_000008.11', 'NC_000009.12', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12', 'NC_000013.11', 'NC_000014.9', 'NC_000015.10', 'NC_000016.10', 'NC_000017.11', 'NC_000018.10', 'NC_000019.10', 'NC_000020.11', 'NC_000021.9', 'NC_000022.11', 'NC_000023.11', 'NC_000024.10']:
                normalizedVariant["b38SPDI"] = spdi["seq_id"] + ":" + str(spdi["position"]) + ":" + spdi["deleted_sequence"] + ":" + spdi["inserted_sequence"]
                normalizedVariant["b38HGVS"] = json.loads(SPDI_hgvs(normalizedVariant["b38SPDI"]).text)["data"]["hgvs"]
    else:
        abort(500, "Unrecognized variant syntax.")
    return normalizedVariant


def translate_terminology(codeSystem, code):
    # validate input parameters
    if codeSystem not in ['http://www.nlm.nih.gov/research/umls/rxnorm', 'http://snomed.info/sct', 'http://hl7.org/fhir/sid/icd-10']:
        return "unrecognized code system"
    response = [{}]

    # process rxnorm input
    if codeSystem == 'http://www.nlm.nih.gov/research/umls/rxnorm':
        MapRxNorm = fetch_concept_map(44872525)
        if MapRxNorm.status_code == 200:
            MapRxNorm = MapRxNorm.json()
            for element in MapRxNorm["group"][0]["element"]:
                if element["code"] == code:
                    response[0]["outcome"] = 'match found'
                    response[0]["code"] = element["target"][0]["code"]
                    response[0]["system"] = 'https://ncithesaurus.nci.nih.gov'
                    response[0]["display"] = element["target"][0]["display"]
                    return response
            response[0]["outcome"] = 'no match found'
            response[0]["system"] = 'https://ncithesaurus.nci.nih.gov'
            return response
        else:
            abort(500, "HAPI server error")

    # process snomed input, return disease ontology code AND medgen
    if codeSystem == 'http://snomed.info/sct':
        Mapsnomed = fetch_concept_map(44947014)
        if Mapsnomed.status_code == 200:
            Mapsnomed = Mapsnomed.json()
            for element in Mapsnomed["group"][0]["element"]:
                if element["code"] == code:
                    response[0]["outcome"] = 'match found'
                    response[0]["code"] = element["target"][0]["code"]
                    response[0]["system"] = 'https://disease-ontology.org'
                    response[0]["display"] = element["target"][0]["display"]
                    break
                response[0]["outcome"] = 'no match found'
                response[0]["system"] = 'https://disease-ontology.org'
        else:
            abort(500, "HAPI server error")

        Mapsnomed = fetch_concept_map(44872524)
        if Mapsnomed.status_code == 200:
            Mapsnomed = Mapsnomed.json()
            response.append({})
            for element in Mapsnomed["group"][0]["element"]:
                if element["code"] == code:
                    response[1]["outcome"] = 'match found'
                    response[1]["code"] = element["target"][0]["code"]
                    response[1]["system"] = 'https://www.ncbi.nlm.nih.gov/medgen'
                    response[1]["display"] = element["target"][0]["display"]
                    break
                response[1]["outcome"] = 'no match found'
                response[1]["system"] = 'https://www.ncbi.nlm.nih.gov/medgen'
        else:
            abort(500, "HAPI server error")
        return response

    # process ICD10 input, return disease ontology AND medgen codes
    if codeSystem == 'http://hl7.org/fhir/sid/icd-10':
        code = code.upper().replace(".", "")
        MapICD10 = fetch_concept_map(44872527)
        if MapICD10.status_code == 200:
            MapICD10 = MapICD10.json()
            for element in MapICD10["group"][0]["element"]:
                if element["code"] == code:
                    response[0]["outcome"] = 'match found'
                    response[0]["code"] = element["target"][0]["code"]
                    response[0]["system"] = 'https://disease-ontology.org'
                    response[0]["display"] = element["target"][0]["display"]
                    break
                response[0]["outcome"] = 'no match found'
                response[0]["system"] = 'https://disease-ontology.org'
        else:
            abort(500, "HAPI server error")

        MapICD10 = fetch_concept_map(44872523)
        if MapICD10.status_code == 200:
            MapICD10 = MapICD10.json()
            response.append({})
            for element in MapICD10["group"][0]["element"]:
                if element["code"] == code:
                    response[1]["outcome"] = "match found"
                    response[1]["code"] = element["target"][0]["code"]
                    response[1]["system"] = "https://www.ncbi.nlm.nih.gov/medgen"
                    response[1]["display"] = element["target"][0]["display"]
                    break
                response[1]["outcome"] = "no match found"
                response[1]["system"] = "https://www.ncbi.nlm.nih.gov/medgen"
        else:
            abort(500, "HAPI server error")
        return response


def normalize_variant_hgvs(variant):
    try:
        return input_normalization.normalize(variant)
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        abort(422, 'Failed LiftOver')


def seqfetcher(acc, start, end):
    try:
        return get_ref_seq_subseq(acc, start, end)
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        abort(404, 'Not Found')


def normalize_hla(allele):
    try:
        return {
            allele: {
                "G": hla.redux(allele, "G"),
                "P": hla.redux(allele, "P"),
                "lg": hla.redux(allele, "lg"),
                "lgx": hla.redux(allele, "lgx"),
                "W": hla.redux(allele, "W"),
                "exon": hla.redux(allele, "exon"),
                "U2": hla.redux(allele, "U2"),
                "S": hla.redux(allele, "S")
            }
        }
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        abort(422, 'Failed HLA normalization')
