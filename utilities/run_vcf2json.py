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

    # Step 1: Extract data into variants_data, molecular_output, phase_output
    extractData("vcfData.csv", variants_data, molecular_output, phase_output)

    # Step 2: Separate wildtypes from variants
    wildtypes_data = []
    filtered_variants_data = []
    bed_entries = []

    current_bed_entry = None

    for variant in variants_data:
        gt = variant.get("GT")
        alt = variant.get("ALT")
        ref = variant.get("REF")
        chrom = variant.get("CHROM")
        pos = variant.get("POS")

        # Determine if the record is a wildtype (non-variant)
        if gt in ['0/0', '0|0', '0']:
            wildtypes_data.append(variant)

            # Create or extend a BED entry for contiguous non-variants
            if current_bed_entry and current_bed_entry["CHROM"] == chrom and current_bed_entry["END"] + 1 == pos - 1:
                # Extend the current bed entry
                current_bed_entry["END"] = pos - 1
            else:
                # Start a new bed entry
                if current_bed_entry:
                    bed_entries.append(current_bed_entry)
                current_bed_entry = OrderedDict({
                    "CHROM": chrom,
                    "START": pos - 1,
                    "END": pos - 1
                })
        else:
            filtered_variants_data.append(variant)
            # Close the current bed entry for non-variants if one is open
            if current_bed_entry:
                bed_entries.append(current_bed_entry)
                current_bed_entry = None

    # Append the last BED entry if it exists
    if current_bed_entry:
        bed_entries.append(current_bed_entry)

    # Step 3: Filter wildtypes_data to include only required fields
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
            "GT": wt.get("GT"),
            "SPDI": wt.get("SPDI"),
            "allelicState": wt.get("allelicState"),
        }
        wildtypes_filtered.append(wt_filtered)

    # Replace wildtypes_data with the filtered data
    wildtypes_data = wildtypes_filtered

    # Step 4: Write the output files
    files_to_write = [
        ("variantsData.json", filtered_variants_data),
        ("molecularConsequences.json", molecular_output),
        ("phaseData.json", phase_output),
        ("wildtypes.json", wildtypes_data),
        ("BED.json", bed_entries)
    ]

    for filename, data in tqdm(files_to_write, desc="Writing files"):
        with open(filename, "w") as f:
            f.write(json.dumps(data, indent=4))

    # Step 5: Remove temporary files if needed
    try:
        os.remove("convertedVCF.json")
    except FileNotFoundError:
        pass  # File doesn't exist, no action needed

    print("Data Generated.")

run_vcf2json()
