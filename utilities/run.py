import vcf2json
import json
import os
import common
import pandas as pd
from tqdm import tqdm
from collections import OrderedDict


def exrtractData(csv_file_path, variants_data, molecular_output, phase_output):
    try:
        # Fetching all Transcripts data and storing it in transcript_map
        transcript_map = {}
        common.query_genes(transcript_map)

        df = pd.read_csv(csv_file_path)
        for _, row in tqdm(df.iterrows(), total=len(df), desc= "Processing Data"):
            phased_rec_map = {}
            vcf2json.vcf2json(row['vcf_filename'],
                              row['ref_build'],
                              row['patient_id'],
                              row['test_date'],
                              row['test_id'],
                              row['specimen_id'],
                              row['genomic_source_class'],
                              row['ratio_ad_dp'],
                              row['sample_position'], transcript_map, variants_data, molecular_output, phased_rec_map)

            file = pd.read_json("convertedVCF.json", orient = str)

            vcf2json.add_phased_relationship_obv(row['patient_id'],
                                                 row['test_id'],
                                                 row['specimen_id'],
                                                 row['ref_build'], phase_output, file, phased_rec_map)

    except Exception as e:
        print("Error with csv file:", e)

def run_vcf2json():
    variants_data = []
    molecular_output = []
    phase_output = []
    exrtractData("vcfData.csv", variants_data, molecular_output, phase_output)

    files_to_write = [
        ("variantsData.json", variants_data),
        ("molecularConsequences.json", molecular_output),
        ("phaseData.json", phase_output)]

    for filename, data in tqdm(files_to_write, desc="Writing files"):
        with open(filename, "w") as f:
            f.write(json.dumps(data, indent=4))

    # Deletion of convertedVCF.JSON file
    try :
        os.remove("convertedVCF.json")
    except:
        print("Error! can't delete vcf file.")

    print("Data Generated.")

run_vcf2json()
