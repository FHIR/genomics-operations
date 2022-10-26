import vcf
import json
from collections import OrderedDict
from gene_ref_seq import _get_ref_seq_by_chrom
from SPDI_Normalization import get_normalized_spdi
from common import *
import copy
import pandas as pd
import re

phased_rec_map = {}
spr_json = []

df = pd.read_json("variants.json", lines=True)

def add_phase_records(record):
    if(record.samples[0].phased is False):
        return
    sample_data = record.samples[0].data
    if(sample_data.GT is not None and
       len(sample_data.GT.split('|')) >= 2 and
       'PS' in sample_data._fields):
        sample_data_ps = sample_data.PS
        if isinstance(sample_data.PS, list):
            sample_data_ps = sample_data_ps[0]
        phased_rec_map.setdefault(sample_data_ps, []).append(record)


def add_phased_relationship_obv(patientID, test_id, specimen_id, ref_build):
    sequence_rels = get_sequence_relation(phased_rec_map)
    c = len(sequence_rels)
    df_func = df[(df['patientID']==patientID)]
    for index in sequence_rels.index:
        output_json = OrderedDict()
        output_json['patientID'] = patientID

        relation = sequence_rels.at[index, 'Relation']
        pos1 = sequence_rels.at[index, 'POS1']
        pos2 = sequence_rels.at[index, 'POS2']

        df_copy = df_func[((df_func['testID']==test_id) & (df_func['specimenID']==specimen_id) & (df_func['genomicBuild']==ref_build)) & ((df_func['POS']==pos1-1) | (df_func['POS']==pos2-1))]

        if len(df_copy) != 2:
            c-=1
            continue

        output_json['variantID1'] = str(df_copy.iloc[0]['_id']['$oid'])
        output_json['variantID2'] = str(df_copy.iloc[1]['_id']['$oid'])
        output_json['phase'] = relation
        
        spr_json.append(output_json)
        c-=1

def _valid_record(record, genomic_source_class, sample_position):
    svAltRegex = re.compile("^<{1}.*>{1}$");
    if len(record.samples) < 1:
        return False
    if not (validate_chrom_identifier(record.CHROM)):
        return False
    if not hasattr(record.samples[sample_position].data, "GT"):
        return False
    if record.is_sv:
        if len(record.samples) > 1:
            return False
        if(record.INFO['SVTYPE'][0].upper() not in list(SVs)):
            return False
        if(not all(alt is None or alt.type in ['SNV', 'MNV'] or
           isinstance(alt, vcf.model._SV) or svAltRegex.match(str(alt)) or (str(alt).isalpha() or (alt == '.' and len(record.ALT) == 1))
           for alt in record.ALT)):
            return False
        if(record.INFO['SVTYPE'][0].upper() in list(SVs - {'DUP', 'CNV'}) and
           '.' in record.samples[sample_position]["GT"] and
           genomic_source_class.lower() == Genomic_Source_Class.GERMLINE.value.lower()):
            return False
    else:
        if(not all(alt is None or ((alt.type in ['SNV', 'MNV'] or '*' not in str(alt)) and str(alt).isalpha()) for alt in record.ALT)):
            return False
        if('.' in record.samples[sample_position]["GT"] and
           genomic_source_class.lower() == Genomic_Source_Class.GERMLINE.value.lower()):
            return False
    if(record.FILTER is not None and len(record.FILTER) != 0):
        return False
    if record.samples[sample_position]["GT"] in ['0/0', '0|0', '0']:
        return False
    if not record.REF.isalpha():
        return False
    if record.CHROM == "M" and (
        (len(
            record.samples[sample_position].gt_alleles) == 1 and
            record.samples[sample_position].gt_alleles[0] == "0") or len(
            record.samples[sample_position].gt_alleles) == 2):
        return False
    return True


def vcf2json(vcf_filename=None, ref_build=None, patient_id=None,
             test_date=None, test_id=None, specimen_id=None,
             genomic_source_class=None, ratio_ad_dp=0.99, sample_position=0):
    
    output_json_array = []
    if not (vcf_filename):
        raise Exception('You must provide vcf_filename')
    if not ref_build or ref_build not in ["GRCh37", "GRCh38"]:
        raise Exception(
            'You must provide build number ("GRCh37" or "GRCh38")')
    if not (patient_id):
        raise Exception('You must provide patient_id')
    if not (test_date):
        raise Exception('You must provide test_date')
    if not (test_id):
        raise Exception('You must provide test_id')
    if not (specimen_id):
        raise Exception('You must provide specimen_id')
    if genomic_source_class is not None and genomic_source_class.title() not in Genomic_Source_Class.set_():
        raise Exception(
            ("Please provide a valid Genomic Source Class " +
             "('germline' or 'somatic' or 'mixed')"))

    try:
        vcf_reader = vcf.Reader(filename=vcf_filename)
    except FileNotFoundError:
        raise
    except BaseException:
        raise Exception("Please provide valid  'vcf_filename'")
    if not patient_id:
        patient_id = vcf_reader.samples[sample_position]

    for record in vcf_reader:
        if not _valid_record(record, genomic_source_class, sample_position):
            continue
        add_phase_records(record)
        output_json = OrderedDict()
        output_json["patientID"] = patient_id
        output_json["testDate"] = test_date
        output_json["testID"] = test_id
        output_json["specimenID"] = specimen_id
        output_json["genomicBuild"] = ref_build
        record.CHROM = extract_chrom_identifier(record.CHROM)
        output_json["CHROM"] = f"chr{record.CHROM}"
        output_json["POS"] = record.POS - 1
        output_json["REF"] = record.REF
        
        alts = record.ALT
        noRefFlag = 0

        if record.FILTER is None:
            output_json["FILTER"] = '.'
        elif isinstance(record.FILTER, list) and len(record.FILTER)==0:
            output_json["FILTER"] = 'PASS'

        if 'SVTYPE' in record.INFO and record.INFO['SVTYPE'] is not None:
            output_json["SVTYPE"] = record.INFO['SVTYPE']
            if record.INFO['SVTYPE'] == 'INS':
                output_json["POS"] = record.POS
        if 'CIPOS' in record.INFO and record.INFO['CIPOS'] is not None:
            output_json["CIPOS"] = record.INFO['CIPOS']
        if 'CIEND' in record.INFO and record.INFO['CIEND'] is not None:
            output_json["CIEND"] = record.INFO['CIEND']
        if 'END' in record.INFO and record.INFO['END'] is not None:
            output_json["END"] = record.INFO['END']


        hasAD = False
        alt_ad_index = 1
        output_json["GT"] = record.samples[sample_position]["GT"]
        if hasattr(record.samples[sample_position].data, "PS") and record.samples[sample_position]["PS"] is not None:
            output_json["PS"] = record.samples[sample_position]["PS"]
        if hasattr(record.samples[sample_position].data, "CN") and record.samples[sample_position]["CN"] is not None:
            output_json["CN"] = record.samples[sample_position]["CN"]
        if hasattr(record.samples[sample_position].data, "AD") and record.samples[sample_position]["AD"] is not None:
            hasAD = True
            output_json["ADS"] = []

            for index in range(0, len(record.samples[sample_position]["AD"])):
                if record.samples[sample_position]["AD"][index] == None:
                    record.samples[sample_position]["AD"][index] = 0

            #Split the genotype into a list of integers. This will help find the number of alternate alleles in the VCF row
            genotypeList = re.split(r'\D+', record.samples[sample_position]["GT"])      

            for index in range(0, len(genotypeList)):
                if genotypeList[index] == (None or ''):
                    genotypeList[index] = 0
                else:
                    genotypeList[index] = int(genotypeList[index])

            #altNumber is the number of alternate alleles in the VCF row
            altNumber = max(genotypeList)

            #VCF row contains refAD count followed by AD counts for each alt allele
            if len(record.samples[sample_position]["AD"]) == altNumber + 1:
                output_json["ADS"].append({"AD": int(record.samples[sample_position]["AD"][0])})
            
            #VCF row contains only alt allele AD counts
            else:
                output_json["ADS"].append({"AD": 0})
                noRefFlag = 1
                #The alt allele ADs start at 0, since there is no reference AD count taking that slot in the VCF row
                alt_ad_index = 0


        if hasattr(record.samples[sample_position].data, "DP") and record.samples[sample_position]["DP"] is not None:
            output_json["DP"] = int(record.samples[sample_position]["DP"])


        ref_seq = _get_ref_seq_by_chrom(ref_build, extract_chrom_identifier(record.CHROM))

        if not record.is_sv and record.ALT is not None and all(alt is not None for alt in record.ALT):
            spdi = (f'{ref_seq}:{record.POS - 1}:{record.REF}:' +
                    f'{"".join(list(map(str, list(record.ALT))))}')
            # Calculate SPDI directly for SNVs and MNVs
            if(all(alt is not None and alt.type in ['SNV', 'MNV'] for alt in record.ALT)):
                output_json["SPDI"] = spdi

            # Calculate SPDI using NCBI API for InDel
            if len(record.REF) != len("".join(list(map(str, list(record.ALT))))):
                output_json["SPDI"] = get_normalized_spdi(ref_seq, (record.POS - 1), record.REF, "".join(list(map(str, list(record.ALT)))), ref_build)            

        alleles = get_allelic_state(record, ratio_ad_dp)

        if (alleles['CODE'] != "" or alleles['ALLELE'] != "") and genomic_source_class.lower() == Genomic_Source_Class.GERMLINE.value.lower() :
            output_json["allelicState"] = alleles['ALLELE']

        output_json["genomicSourceClass"] = genomic_source_class
        onRef = 1

        extractINFOField(output_json, record, codeDict)

        for alt in alts:
            altDict, alt_ad_index, onRef = getMultADs(output_json, record, sample_position, alt, alt_ad_index, hasAD, noRefFlag, onRef)
            output_json_array.append(altDict)

    output_json_string = json.dumps(output_json_array, indent=4)

    fileOutput = open("convertedVCF.json", "w")
    fileOutput.write(output_json_string)
    fileOutput.close()

    return output_json_array

#getMultADs finds the allele reads for each decomposed alternate allele in the VCF row. 
def getMultADs(output_json, record, sample_position, alt, alt_ad_index, hasAD, noRefFlag, onRef):
    altDict = OrderedDict()
    if alt:
        if "ALT" in output_json:
            del output_json["ALT"]

        for key, value in output_json.items():
            altDict[key]=copy.deepcopy(value)
            if key == 'REF':
                altDict["ALT"] = f"{alt}"

        try:
            if hasAD and noRefFlag == 0:
                if record.samples[sample_position]["AD"][alt_ad_index] != '.' and record.samples[sample_position]["AD"][alt_ad_index] is not None:
                    altDict["ADS"].append({"AD": int(record.samples[sample_position]["AD"][alt_ad_index])})
                else:
                    altDict["ADS"].append({"AD": 0})
                alt_ad_index += 1
            elif noRefFlag == 1:
                #If there is no reference AD, set initial AD to 0 and insert new ADs instead of appending
                if record.samples[sample_position]["AD"][alt_ad_index] != '.' and record.samples[sample_position]["AD"][alt_ad_index] is not None:
                    altDict["ADS"].insert(0, {"AD": int(record.samples[sample_position]["AD"][alt_ad_index])})
                else:
                    altDict["ADS"].insert(0, {"AD": 0})
                alt_ad_index += 1
        except Exception as e:
            altDict["ADS"].append({"AD": 0})
            alt_ad_index += 1
    else:
        if "ALT" in output_json:
            del output_json["ALT"]
        for key, value in output_json.items():
            altDict[key]=copy.deepcopy(value)


    return altDict, alt_ad_index, onRef

#Processes SnpEff and SnpSift output into final json file.
def extractINFOField(output_json, record, codeDict):
    #Info field is dict with entries: ANN (SnpEff output) and POPAF (gnomAD output)
    
    if 'ANN' in record.INFO:
        output_json["molecConseq"] = []
        firstFlag = True
        for ann in record.INFO['ANN']:
            parseANN(output_json, ann, firstFlag, codeDict)
            firstFlag = False

    if 'POPAF' in record.INFO:
        for popAF in record.INFO['POPAF']:
            if popAF is not None:
                output_json["popAlleleFreq"] = float(popAF)

    if 'LOF' in record.INFO:
        output_json["funcConseq"] = []
        output_json["funcConseq"].append({"system": r'http://sequenceontology.org/',
                                            "code": "SO_0002054",
                                            "display": "loss_of_function_variant"})

#Orders and extracts molecular consequence data from SnpEff annotations.
def parseANN(output_json, ann, firstFlag, codeDict):
    annList = ann.split('|')
    conseqList = annList[1].split('&')

    for conseq in conseqList:
        uniqueConseq = True
        for conseqIterator in output_json["molecConseq"]:
            if conseqIterator["display"] == conseq:
                uniqueConseq = False
        
        if uniqueConseq:
            try:
                system = codeDict[conseq][0]
                code = codeDict[conseq][1]
            
                output_json["molecConseq"].append({"system": system,
                                                    "code": code,
                                                    "display": conseq})
            except:
                print("molecular consequence: " + conseq + " not represented in code table")
 
    if firstFlag:    
        output_json["predictedMolecImpact"] = annList[2]
