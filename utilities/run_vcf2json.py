import vcf2json
import json
import os
import common
import pandas as pd
from tqdm import tqdm

def extractData(csv_file_path, variants_data, molecular_output, phase_output):
    try:
        # Fetching all Transcripts data and storing it in transcript_map
        transcript_map = {}
        common.query_genes(transcript_map)

        df = pd.read_csv(csv_file_path)

        for _, row in tqdm(df.iterrows(), total=len(df)-1, desc="Processing Data"):
            phased_rec_map = {}
            vcf2json.vcf2json(
                row['vcf_filename'],
                row['ref_build'],
                row['patient_id'],
                row['test_date'],
                row['test_id'],
                row['specimen_id'],
                row['genomic_source_class'],
                row['ratio_ad_dp'],
                row['sample_position'],
                transcript_map,
                variants_data,
                molecular_output,
                phased_rec_map
            )

            file = pd.read_json("convertedVCF.json", orient=str)

            vcf2json.add_phased_relationship_obv(
                row['patient_id'],
                row['test_id'],
                row['specimen_id'],
                row['ref_build'],
                phase_output,
                file,
                phased_rec_map
            )

    except Exception as e:
        print("Error with csv file:", e)

def run_vcf2json():
    variants_data = []
    molecular_output = []
    phase_output = []

    extractData("vcfData.csv", variants_data, molecular_output, phase_output)

    # Separate wildtypes from variants
    wildtypes_data = []
    filtered_variants_data = []

    for variant in variants_data:
        gt = variant.get("GT")
        alt = variant.get("ALT")
        ref = variant.get("REF")

        # Determine if the record is a wildtype (non-variant)
        if gt in ['0/0', '0|0', '0'] or alt == ref or alt in [None, '.', '']:
            wildtypes_data.append(variant)
        else:
            filtered_variants_data.append(variant)

    # Update variants_data to contain only variants
    variants_data = filtered_variants_data

    # Filter wildtypes_data to include only required fields
    wildtypes_filtered = []
    for wt in wildtypes_data:
        # Create a new dictionary with only the required fields
        wt_filtered = {
            "_id": wt.get("_id"),  # Unique identifier
            "patientID": wt.get("patientID"),
            "testDate": wt.get("testDate"),
            "testID": wt.get("testID"),
            "specimenID": wt.get("specimenID"),
            "genomicBuild": wt.get("genomicBuild"),
            "CHROM": wt.get("CHROM"),
            "POS": wt.get("POS"),
            "REF": wt.get("REF"),
            "ALT": wt.get("ALT"),  # May be null for wildtypes
            "END": wt.get("END"),  # Optional
            "FILTER": wt.get("FILTER"),
            "GT": wt.get("GT"),  # Optional
            "SPDI": wt.get("SPDI"),  # Optional
            "allelicState": wt.get("allelicState"),  # Optional
            "geneCode": wt.get("geneCode"),  # Optional
            "geneDesc": wt.get("geneDesc"),  # Optional
        }
        wildtypes_filtered.append(wt_filtered)

    # Replace wildtypes_data with the filtered data
    wildtypes_data = wildtypes_filtered

    # Proceed to write the output files
    files_to_write = [
        ("variantsData.json", variants_data),
        ("molecularConsequences.json", molecular_output),
        ("phaseData.json", phase_output),
        ("wildtypes.json", wildtypes_data)
    ]

    for filename, data in tqdm(files_to_write, desc="Writing files"):
        with open(filename, "w") as f:
            f.write(json.dumps(data, indent=4))

    # Remove temporary files if needed
    try:
        os.remove("convertedVCF.json")
    except FileNotFoundError:
        pass  # File doesn't exist, no action needed

    print("Data Generated.")

run_vcf2json()
