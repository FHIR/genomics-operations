import vcf2json
import json
import os
import common
import pandas as pd
from tqdm import tqdm
from collections import OrderedDict
import parse_non_variant

def extractData(csv_file_path, variants_data, molecular_output, phase_output):
    try:
        # Fetching all Transcripts data and storing it in transcript_map
        transcript_map = {}
        common.query_genes(transcript_map)

        df = pd.read_csv(csv_file_path)
        for _, row in tqdm(df.iterrows(), total=len(df)-1, desc="Processing Data"):
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

            file = pd.read_json("convertedVCF.json", orient=str)

            vcf2json.add_phased_relationship_obv(row['patient_id'],
                                                 row['test_id'],
                                                 row['specimen_id'],
                                                 row['ref_build'], phase_output, file, phased_rec_map)

    except Exception as e:
        print("Error with csv file:", e)


def extract_non_variant_data(csv_file_path):
    """
    Extracts non-variant data from a CSV file and generates combined JSON outputs.

    Args:
        csv_file_path (str): Path to the CSV file containing GVCF metadata.

    Returns:
        tuple: Lists containing non-variant JSON data and BED JSON entries.
    """
    df = pd.read_csv(csv_file_path)
    processed_files = set()  # To track processed GVCF files
    bed_entries = []  # Collect BED entries for all files
    non_variant_entries = []  # Collect non-variant entries for all files

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing GVCF Data"):
        gvcf_filename = row['vcf_filename']
        ref_build = row['ref_build']  # Dynamically set build version

        # Skip already processed files
        if gvcf_filename in processed_files:
            continue

        # Extract BED entries and detailed records from the GVCF file
        if gvcf_filename.endswith(".g.vcf"):
            try:
                studied_regions, detailed_blocks = parse_non_variant.parse_gvcf(
                    gvcf_filename, int(ref_build[-2:])
                )  # GRCh37 or GRCh38
            except Exception as e:
                print(f"Error processing {gvcf_filename}: {e}")
                continue

            # Update studied regions (BED data) with metadata
            for region in studied_regions:
                region.update({
                    "REF_BUILD": ref_build,
                    "PATIENT_ID": row['patient_id'],
                    "TEST_DATE": row['test_date'],
                    "TEST_ID": row['test_id'],
                    "SPECIMEN_ID": row['specimen_id'],
                    "GENOMIC_SOURCE_CLASS": row['genomic_source_class'],
                    "RATIO_AD_DP": row['ratio_ad_dp'],
                    "SAMPLE_POSITION": row['sample_position'],
                })
                bed_entries.append(region)

            # Update detailed blocks with metadata
            for block in detailed_blocks:
                block.update({
                    "REF_BUILD": ref_build,
                    "PATIENT_ID": row['patient_id'],
                    "TEST_DATE": row['test_date'],
                    "TEST_ID": row['test_id'],
                    "SPECIMEN_ID": row['specimen_id'],
                    "GENOMIC_SOURCE_CLASS": row['genomic_source_class'],
                    "RATIO_AD_DP": row['ratio_ad_dp'],
                    "SAMPLE_POSITION": row['sample_position'],
                })
                non_variant_entries.append(block)

            # Mark file as processed
            processed_files.add(gvcf_filename)

    return non_variant_entries, bed_entries


def transform_non_variant_output(non_variant_data, bed_file_output):
    """
    Transforms non-variant data to create a record for each position, keeping only the single allele at each position.

    Args:
        non_variant_data (list): List of non-variant entries with START and REF_ALLELE.

    Returns:
        list: Transformed list with a record for each position and its single allele.
    """
    transformed_data = []
    transformed_data2 = []

    for record in tqdm(non_variant_data, desc="Transforming Non-Variant Data"):
        chrom = record["CHROM"]
        ref_allele = record["REF_ALLELE"]
        FILTER = record["FILTER"]
        ref_build = record["REF_BUILD"]
        patient_id = record["PATIENT_ID"]
        test_date = record["TEST_DATE"]
        test_id = record["TEST_ID"]
        specimen_id = record["SPECIMEN_ID"]
        genomic_source_class = record["GENOMIC_SOURCE_CLASS"]
        ratio_ad_dp = record["RATIO_AD_DP"]
        gt = record["GT"]
        allelicState = record["allelicState"]
        pos = record["POS"]  # Take the starting position directly

        # Create a record for each position using REF_ALLELE length
        for idx, allele in enumerate(ref_allele):
            transformed_data.append({
                "CHROM": chrom,
                "REF": allele,# Single allele
                "ALT": allele,
                "POS": pos + idx,  # Increment position
                "GT": gt,
                "FILTER": FILTER,
                "allelicState": allelicState,
                "genomicBuild": ref_build,
                "patientID": patient_id,
                "testDate": test_date,
                "testID": test_id,
                "specimenID": specimen_id,
            })

        # Process BED file data (directly use START and END)
    for record in tqdm(bed_file_output, desc="Transforming BED File Data"):
        chrom = record["CHROM"]
        start = record["START"]
        end = record["END"]
        patient_id = record["PATIENT_ID"]

        transformed_data2.append({
            "Chrom": chrom,
            "Start": start,
            "End": end,
            "patientID": patient_id,
        })

    return transformed_data, transformed_data2


def run_vcf2json():
    variants_data = []
    molecular_output = []
    phase_output = []

    # Set to True to include non-variant JSON output
    allow_non_variants = True
    files_to_write = []

    if allow_non_variants:
        non_variant_output, bed_file_output = extract_non_variant_data("vcfData.csv")

        # Transform non-variant data into the desired format
        transformed_data, transformed_data2 = transform_non_variant_output(non_variant_output, bed_file_output)

        # Add non-variant files to write list
        files_to_write.extend([
            ("non_variant_output.json", transformed_data),  # Use transformed data
            ("bed_file_output.json", transformed_data2)
        ])

    extractData("vcfData.csv", variants_data, molecular_output, phase_output)

    # Add variant data files to the write list
    files_to_write.extend([
        ("variantsData.json", variants_data),
        ("molecularConsequences.json", molecular_output),
        ("phaseData.json", phase_output),
    ])

    # Write files
    for filename, data in tqdm(files_to_write, desc="Writing files"):
        with open(filename, "w") as f:
            f.write(json.dumps(data, indent=4))

    # Clean up temporary files
    try:
        os.remove("convertedVCF.json")
    except FileNotFoundError:
        pass

    print("Data Generated.")


if __name__ == "__main__":
    run_vcf2json()
