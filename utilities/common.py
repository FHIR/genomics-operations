import os
import re
from enum import Enum

import pandas as pd
import pymongo

from pathlib import Path
from dotenv import load_dotenv

BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent
load_dotenv(PROJECT_ROOT / "secrets.env")

utilities_data_client_uri = f"mongodb+srv://readonly:{os.getenv('MONGODB_READONLY_PASSWORD')}@cluster0.8ianr.mongodb.net/UtilitiesData"
utilities_client = pymongo.MongoClient(utilities_data_client_uri)
utilities_db = utilities_client.UtilitiesData
transcript_data = utilities_db.Transcripts

GERMLINE = 'Germline'
SOMATIC = 'Somatic'
MIXED = 'Mixed'
SVs = {'INS', 'DEL', 'DUP', 'CNV', 'INV'}

# Dictionary of SO Codes for use in molecular consequence pipeline
codeDict = {'coding_sequence_variant ': ['http://sequenceontology.org', 'SO:0001580'],
            'chromosome': ['http://sequenceontology.org', 'SO:0000340'],
            'duplication': ['http://sequenceontology.org', 'SO:1000035'],
            'inversion': ['http://sequenceontology.org', 'SO:1000036'], 'inframe_insertion': ['http://sequenceontology.org', 'SO:0001821'], 'disruptive_inframe_insertion': ['http://sequenceontology.org', 'SO:0001824'], 'inframe_deletion': ['http://sequenceontology.org', 'SO:0001822'], 'disruptive_inframe_deletion': ['http://sequenceontology.org', 'SO:0001826'], 'downstream_gene_variant': ['http://sequenceontology.org', 'SO:0001632'], 'exon_variant': ['http://sequenceontology.org', 'SO:0001791'], 'exon_loss_variant': ['http://sequenceontology.org', 'SO:0001572'], 'frameshift_variant': ['http://sequenceontology.org', 'SO:0001589'], 'gene_variant': ['http://sequenceontology.org', 'SO:0001564'], 'feature_ablation': ['http://sequenceontology.org', 'SO:0001879'], 'gene_fusion': ['http://sequenceontology.org', 'SO:0001565'], 'bidirectional_gene_fusion': ['http://sequenceontology.org', 'SO:0002086'], 'rearranged_at_DNA_level': ['http://sequenceontology.org', 'SO:0000904'], 'intergenic_region': ['http://sequenceontology.org', 'SO:0000605'], 'conserved_intergenic_variant': ['http://sequenceontology.org', 'SO:0002017'], 'intragenic_variant': ['http://sequenceontology.org', 'SO:0002011'], 'intron_variant': ['http://sequenceontology.org', 'SO:0001627'], 'conserved_intron_variant': ['http://sequenceontology.org', 'SO:0002018'], 'miRNA': ['http://sequenceontology.org', 'SO:0000276'], 'missense_variant': ['http://sequenceontology.org', 'SO:0001583'], 'initiator_codon_variant': ['http://sequenceontology.org', 'SO:0001582'], 'stop_retained_variant': ['http://sequenceontology.org', 'SO:0001567'], 'protein_protein_contact': ['http://sequenceontology.org', 'SO:0001093'], 'structural_interaction_variant': ['http://sequenceontology.org', 'SO:0002093'], 'rare_amino_acid_variant': ['http://sequenceontology.org', 'SO:0002008'], 'splice_acceptor_variant': ['http://sequenceontology.org', 'SO:0001574'], 'splice_donor_variant': ['http://sequenceontology.org', 'SO:0001575'], 'splice_region_variant': ['http://sequenceontology.org', 'SO:0001630'], 'stop_lost': ['http://sequenceontology.org', 'SO:0001578'], '5_prime_UTR_premature_start_codon_gain_variant': ['http://sequenceontology.org', 'SO:0001988'], 'start_lost': ['http://sequenceontology.org', 'SO:0002012'], 'stop_gained': ['http://sequenceontology.org', 'SO:0001587'], 'synonymous_variant': ['http://sequenceontology.org', 'SO:0001819'], 'start_retained': ['http://sequenceontology.org', 'SO:0002019'], 'transcript_variant': ['http://sequenceontology.org', 'SO:0001576'], 'regulatory_region_variant': ['http://sequenceontology.org', 'SO:0001566'], 'upstream_gene_variant': ['http://sequenceontology.org', 'SO:0001631'], '3_prime_UTR_variant': ['http://sequenceontology.org', 'SO:0001624'], '3_prime_UTR_truncation': ['http://sequenceontology.org', 'SO:0002015'], 'exon_loss': ['http://sequenceontology.org', 'SO:0001572'], '5_prime_UTR_variant': ['http://sequenceontology.org', 'SO:0001623'], '5_prime_UTR_truncation': ['http://sequenceontology.org', 'SO:0002013'], 'sequence_feature': ['http://sequenceontology.org', 'SO:0000110'], 'non_coding_transcript_exon_variant': ['http://sequenceontology.org', 'SO:0001792'], 'conservative_inframe_insertion': ['http://sequenceontology.org', 'SO:0001823'], 'non_coding_transcript_variant': ['http://sequenceontology.org', 'SO:0001619'], 'conservative_inframe_deletion': ['http://sequenceontology.org', 'SO:0001825'],
            'start_retained_variant': ['http://sequenceontology.org', 'SO:0002019']}


class Genomic_Source_Class(Enum):

    @classmethod
    def set_(cls):
        return set(map(lambda c: c.value, cls))

    GERMLINE = GERMLINE
    SOMATIC = SOMATIC
    MIXED = MIXED


def validate_chrom_identifier(chrom):
    chrom = extract_chrom_identifier(chrom)
    pattern = '^[1-9]$|^1[0-9]$|^2[0-2]$|^[XYM]$'
    result = re.match(pattern, chrom)
    return bool(result)


def _get_chrom(chrom_index):
    switcher = {
        23: 'X',
        24: 'Y',
        25: 'M'
    }
    return switcher.get(chrom_index, str(chrom_index))


def extract_chrom_identifier(chrom):
    chrom = chrom.upper().replace("CHR", "")
    if chrom == "MT":
        chrom = "M"
    return chrom


def get_allelic_state(record, ratio_ad_dp, sample_position):
    allelic_state = ''
    allelic_code = ''
    try:
        allelic_frequency = float(record.aaf[0])
    except ZeroDivisionError:
        allelic_frequency = 0
    sample = record.samples[sample_position]
    alleles = sample.gt_alleles
    if record.CHROM != 'M':
        # gt_type: hom_ref = 0; het = 1; hom_alt = 2; uncalled = None
        if len(alleles) >= 2 and sample.gt_type == 1:
            allelic_state = 'heterozygous'
            allelic_code = 'LA6706-1'
        elif len(alleles) >= 2:
            allelic_state = 'homozygous'
            allelic_code = 'LA6705-3'
        elif sample.gt_type is not None and len(alleles) == 1:
            allelic_state = 'hemizygous'
            allelic_code = 'LA6707-9'
        else:
            _error_log_allelicstate(record)
    elif (sample.gt_type is not None and
          len(alleles) == 1 and
          alleles[0] == '1'):
        if allelic_frequency > ratio_ad_dp:
            allelic_state = "homoplasmic"
            allelic_code = "LA6704-6"
        else:
            allelic_state = "heteroplasmic"
            allelic_code = "LA6703-8"
    else:
        _error_log_allelicstate(record)
    return {
        'ALLELE': allelic_state,
        'CODE': allelic_code,
        'FREQUENCY': allelic_frequency
    }


# got this function from https://github.com/elimuinformatics/vcf2fhir/blob/master/vcf2fhir/common.py
def get_sequence_relation(phased_rec_map):
    Relation_table = pd.DataFrame(columns=['POS1', 'POS2', 'Relation'])
    for key in phased_rec_map:
        prev_record = None
        for record in phased_rec_map[key]:
            if prev_record is None:
                prev_record = record
                continue
            prev_data = prev_record.samples[0].data
            record_data = record.samples[0].data
            if (prev_data.PS == record_data.PS):
                if prev_data.GT == record_data.GT:
                    Relation_table = Relation_table.append(
                        {
                            'POS1': prev_record.POS,
                            'POS2': record.POS,
                            'Relation': 'Cis'
                        },
                        ignore_index=True
                    )
                else:
                    Relation_table = Relation_table.append(
                        {
                            'POS1': prev_record.POS,
                            'POS2': record.POS,
                            'Relation': 'Trans'
                        },
                        ignore_index=True
                    )
            prev_record = record
    return Relation_table


def query_genes(transcript_map):

    query_string = [{'$match': {}},  # Match all documents
                    {'$project': {'_id': 0, 'ncbiGeneSymbol': 0, 'featureType': 0, 'build37RefSeq': 0, 'build37Start': 0,
                                  'build37End': 0, 'build38RefSeq': 0, 'build38Start': 0, 'build38End': 0}}]

    result = list(transcript_data.aggregate(query_string))
    for res in result:
        transcript_map[res['transcriptRefSeq']] = res['MANE']


def _error_log_allelicstate(record):
    pass
