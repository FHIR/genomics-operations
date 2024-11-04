import os
import subprocess
import pysam
import json

# Define S3 path and local download path for the FASTQ file
s3_bucket = "s3://1000genomes/"
fastq_path = "phase3/data/HG00096/sequence_read/SRR062634.filt.fastq.gz"  # Use the actual path
local_fastq = "SRR062634.filt.fastq.gz"
output_json = "output.json"  # Output JSON file

# Step 1: Download the FASTQ file from S3
def download_fastq_from_s3():
    try:
        print("Downloading FASTQ file from S3...")
        subprocess.run(["aws", "s3", "cp", s3_bucket + fastq_path, local_fastq], check=True)
        print("Download completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading file from S3: {e}")

# Step 2: Convert FASTQ file to JSON format
def fastq_to_json(fastq_file, json_file):
    entries = []  # List to hold each entry as a dictionary
    try:
        with pysam.FastxFile(fastq_file) as fastq:
            for entry in fastq:
                # Create a dictionary for each entry with relevant details
                entry_dict = {
                    "name": entry.name,         # Read identifier
                    "sequence": entry.sequence, # Sequence data
                    "quality": entry.quality    # Quality scores (only for FASTQ files)
                }
                entries.append(entry_dict)

        # Write the list of entries to the JSON file
        with open(json_file, "w") as json_out:
            json.dump(entries, json_out, indent=4)

        print(f"JSON file created: {json_file}")
    except Exception as e:
        print(f"Error processing FASTQ to JSON: {e}")

# Step 3: Execute the steps
def main():
    download_fastq_from_s3()
    fastq_to_json(local_fastq, output_json)

    # Clean up by deleting the downloaded FASTQ file
    try:
        os.remove(local_fastq)
        print("Downloaded FASTQ file deleted.")
    except OSError as e:
        print(f"Error deleting FASTQ file: {e}")

if __name__ == "__main__":
    main()
