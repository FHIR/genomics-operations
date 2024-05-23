import vcf
import json
from collections import OrderedDict
from gene_ref_seq import _get_ref_seq_by_chrom
from SPDI_Normalization import get_normalized_spdi
import common
import re
import uuid


def add_phase_records(record, phased_rec_map, sample_position):
    if (record.samples[sample_position].phased is False):
        return
    sample_data = record.samples[sample_position].data
    if (sample_data.GT is not None and
       len(sample_data.GT.split('|')) == 2 and
       sample_data.GT.split('|')[0] != sample_data.GT.split('|')[1] and
       'PS' in sample_data._fields):
        sample_data_ps = sample_data.PS
        if isinstance(sample_data.PS, list):
            sample_data_ps = sample_data_ps[0]
        phased_rec_map.setdefault(sample_data_ps, []).append(record)


def add_phased_relationship_obv(patientID, test_id, specimen_id, ref_build, phase_data, df, phased_rec_map):
    sequence_rels = common.get_sequence_relation(phased_rec_map)
    c = len(sequence_rels)
    df_func = df[(df['patientID'] == patientID)]
    for index in sequence_rels.index:
        output_json = OrderedDict()
        output_json['patientID'] = patientID

        relation = sequence_rels.at[index, 'Relation']
        pos1 = sequence_rels.at[index, 'POS1']
        pos2 = sequence_rels.at[index, 'POS2']

        df_copy = df_func[((df_func['testID'] == test_id) & (df_func['specimenID'] == specimen_id) & (
            df_func['genomicBuild'] == ref_build)) & ((df_func['POS'] == pos1-1) | (df_func['POS'] == pos2-1))]
        if len(df_copy) != 2:
            c -= 1
            continue

        output_json['variantID1'] = str(df_copy.iloc[0]['_id'])
        output_json['variantID2'] = str(df_copy.iloc[1]['_id'])
        output_json['phase'] = relation
        phase_data.append(output_json)
        c -= 1


def _valid_record(record, genomic_source_class, sample_position):
    svAltRegex = re.compile("^<{1}.*>{1}$")
    if len(record.samples) < 1:
        return False
    if not (common.validate_chrom_identifier(record.CHROM)):
        return False
    if not hasattr(record.samples[sample_position].data, "GT"):
        return False
    if record.is_sv:
        if len(record.samples) > 1:
            return False
        if (record.INFO['SVTYPE'].upper() not in list(common.SVs)):
            return False
        if (not all(alt is None or alt.type in ['SNV', 'MNV'] or
           isinstance(alt, vcf.model._SV) or svAltRegex.match(str(alt)) or (str(alt).isalpha() or (alt == '.' and len(record.ALT) == 1))
           for alt in record.ALT)):
            return False
        if (record.INFO['SVTYPE'].upper() in list(common.SVs - {'DUP', 'CNV'}) and
           '.' in record.samples[sample_position]["GT"] and
                genomic_source_class.lower() == common.Genomic_Source_Class.GERMLINE.value.lower()):
            return False
    else:
        if (not all(alt is None or ((alt.type in ['SNV', 'MNV'] or '*' not in str(alt)) and str(alt).isalpha()) for alt in record.ALT)):
            return False
        if ('.' in record.samples[sample_position]["GT"] and
           genomic_source_class.lower() == common.Genomic_Source_Class.GERMLINE.value.lower()):
            return False
    if (record.FILTER is not None and len(record.FILTER) != 0):
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
             genomic_source_class=None, ratio_ad_dp=0.99, sample_position=0,
             transcript_map=None, variants_data=None, molecular_output=None,
             phased_rec_map=None):

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
    if genomic_source_class is not None and genomic_source_class.title() not in common.Genomic_Source_Class.set_():
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
        # high level logic flow:
        #  - get a VCF row using the pyvcf reader
        #  - validate that it is a record we parse
        #  - update the record if it contains multiple ALT alleles
        #  - for each ALT allele in a record, output_json
        #  - write output_json (which is used to compute phase data)
        if not _valid_record(record, genomic_source_class, sample_position):
            continue

        # Convert multi-ALT record to single-ALT record where genotype only includes one of the ALT alleles
        # (example of VCF row that is converted: chr1 236539074 . G A,C 912.77 PASS AC=1 GT:AD:DP 0/2:1,19,9:29)
        # (example of VCF row that is not converted: chr1 236539077 . G A,C 912.77 PASS AC=1 GT:PS:AD:DP 1/2:123:1,19,9:29)
        multiALT = False
        if len(record.ALT) > 1:
            multiALT = True
            if record.samples[sample_position]["GT"] in ["0/1", "0|1", "1/0", "1|0"]:
                record.ALT = [record.ALT[0]]
                GT = record.samples[sample_position]["GT"]
            elif record.samples[sample_position]["GT"] in ["0/2", "0|2", "2/0", "2|0"]:
                record.ALT = [record.ALT[1]]
                GT = record.samples[sample_position]["GT"].replace("2", "1")
            elif record.samples[sample_position]["GT"] in ["1/1", "1|1"]:
                record.ALT = [record.ALT[1]]
            elif record.samples[sample_position]["GT"] in ["2/2", "2|2"]:
                record.ALT = [record.ALT[1]]
                GT = record.samples[sample_position]["GT"].replace("2", "1")

        # Create output record for each remaining ALT allele.
        i = 0
        for alt in record.ALT:
            output_json = OrderedDict()
            if not multiALT:
                add_phase_records(record, phased_rec_map, sample_position)
            output_json["_id"] = uuid.uuid4().hex
            output_json["patientID"] = patient_id
            output_json["testDate"] = test_date
            output_json["testID"] = test_id
            output_json["specimenID"] = specimen_id
            output_json["genomicBuild"] = ref_build
            record.CHROM = common.extract_chrom_identifier(record.CHROM)
            output_json["CHROM"] = f"chr{record.CHROM}"
            output_json["POS"] = record.POS - 1
            output_json["REF"] = record.REF
            # populate ALT
            if len(record.ALT) > 1:
                output_json["ALT"] = str(record.ALT[i])
            else:
                output_json["ALT"] = str(record.ALT[0])

            output_json["END"] = (record.POS - 1 + len(record.REF))
            if record.FILTER is None:
                output_json["FILTER"] = '.'
            elif isinstance(record.FILTER, list) and len(record.FILTER) == 0:
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
            if hasattr(record.samples[sample_position].data, "PS") and record.samples[sample_position]["PS"] is not None:
                output_json["PS"] = record.samples[sample_position]["PS"]
            if hasattr(record.samples[sample_position].data, "CN") and record.samples[sample_position]["CN"] is not None:
                output_json["CN"] = record.samples[sample_position]["CN"]
            # populate GT, allelicFrequency. We are currently setting them to null for multiALT VCF rows
            output_json["GT"] = record.samples[sample_position]["GT"]
            try:
                output_json["allelicFrequency"] = float(record.aaf[0])
            except ZeroDivisionError:
                output_json["allelicFrequency"] = 0
            if multiALT:
                output_json["GT"] = None
                output_json["allelicFrequency"] = 0
                if len(record.ALT) == 1:
                    output_json["GT"] = GT
            output_json["genomicSourceClass"] = genomic_source_class
            # populate SPDI
            ref_seq = _get_ref_seq_by_chrom(ref_build, common.extract_chrom_identifier(record.CHROM))

            if not record.is_sv and output_json["ALT"] is not None:
                spdi = (f'{ref_seq}:{record.POS - 1}:{record.REF}:{output_json["ALT"]}')
                # Calculate SPDI directly for SNVs and MNVs
                if (alt.type in ['SNV', 'MNV']):
                    output_json["SPDI"] = spdi

                # Calculate SPDI using SPDI_Normalization logic
                if len(record.REF) != len(output_json["ALT"]):
                    output_json["SPDI"] = get_normalized_spdi(ref_seq, (record.POS - 1), record.REF, output_json["ALT"], ref_build)
            # populate allelicState
            alleles = common.get_allelic_state(record, ratio_ad_dp, sample_position)

            if (alleles['CODE'] != "" or alleles['ALLELE'] != "") and genomic_source_class.lower() == common.Genomic_Source_Class.GERMLINE.value.lower():
                output_json["allelicState"] = alleles['ALLELE']

            # extractINFOField gathers molecular consequences.  We are currently skipping multiALT VCF rows
            if not multiALT:
                extractINFOField(output_json["_id"], patient_id, record, common.codeDict, molecular_output, output_json, transcript_map)
            variants_data.append(output_json)
            output_json_array.append(output_json)
            i = i + 1

    with open("convertedVCF.json", "w") as f:
        f.write(json.dumps(output_json_array, indent=4))


# Processes SnpEff and SnpSift output into final json file.


def extractINFOField(variant_id, patient_id, record, codeDict, mol_output, output_json, transcript_map):
    # Info field is dict with entries: ANN (SnpEff output) and POPAF (gnomAD output)
    if 'ANN' in record.INFO:
        count = 3
        for ann in record.INFO['ANN']:
            annList = ann.split('|')
            isMane = False

            if annList[6] in transcript_map:
                mane = transcript_map[annList[6]]
                if (mane == 0):
                    isMane = False
                else:
                    isMane = True

            # Checking if MANE or not, If MANE is true for first MolecularConsequences then
            # only one Value is added and others are ignored
            # If MANE is true for 2nd MolecularConsequence add first and second only.
            # If MANE is true for 3rd MolecularConsequences add first, second and third only.
            # If MANE is false for 3rd and true for other then we will add that and break.
            # If no MANE Present only 2 MolecularConsequences are added.

            if (count != 0):
                if (count == 3 and isMane):
                    additionInMolecularConseq(variant_id, patient_id, record, codeDict, mol_output, annList, isMane)
                    break
                elif (count == 2 and isMane):
                    additionInMolecularConseq(variant_id, patient_id, record, codeDict, mol_output, annList, isMane)
                    break
                elif (count == 1 and isMane):
                    additionInMolecularConseq(variant_id, patient_id, record, codeDict, mol_output, annList, isMane)
                    break
                elif (count != 1):
                    additionInMolecularConseq(variant_id, patient_id, record, codeDict, mol_output, annList, isMane)
                    count -= 1

    if 'POPAF' in record.INFO:
        for popAF in record.INFO['POPAF']:
            if popAF is not None:
                output_json["popAlleleFreq"] = float(popAF)

# Orders and extracts molecular consequence data from SnpEff annotations.


def parseANN(molecular_json, annList, firstFlag, codeDict):
    conseqList = annList[1].split('&')
    molecular_json["transcriptRefSeq"] = annList[6]
    molecular_json["cHGVS"] = annList[6] + ":" + annList[9]

    if (annList[10] != ""):
        molecular_json["pHGVS"] = annList[10]

    for conseq in conseqList:
        uniqueConseq = True
        for conseqIterator in molecular_json["featureConsequence"]:
            if conseqIterator["display"] == conseq:
                uniqueConseq = False

        if uniqueConseq:
            try:
                system = codeDict[conseq][0]
                code = codeDict[conseq][1]
                molecular_json["featureConsequence"].append({"system": system,
                                                             "code": code,
                                                             "display": conseq})

            except Exception:
                print("feature consequence: " + conseq + " not represented in code table")

    if firstFlag:
        molecular_json["impact"] = annList[2]


def additionInMolecularConseq(variant_id, patient_id, record, codeDict, mol_output, annList, isMane):
    molecular_json = OrderedDict()
    molecular_json["patientID"] = patient_id
    molecular_json["variantID"] = variant_id

    molecular_json["transcriptRefSeq"] = ""
    molecular_json["MANE"] = isMane
    molecular_json["source"] = "snpEff"
    molecular_json["cHGVS"] = ""
    molecular_json["pHGVS"] = ""
    molecular_json["featureConsequence"] = []

    firstFlag = True

    parseANN(molecular_json, annList, firstFlag, codeDict)

    firstFlag = False

    if 'LOF' in record.INFO:
        molecular_json["functionalEffect"] = []
        molecular_json["functionalEffect"].append({"system": r'http://sequenceontology.org/',
                                                   "code": "SO:0002054",
                                                   "display": "loss_of_function_variant"})

    mol_output.append(molecular_json)
