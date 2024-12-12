import json
import vcf
from tqdm import tqdm
from ref import extract_chromosome_sequence


def normalize_chromosome(chrom):
    """
    Normalizes chromosome names to the `Chr` format.
    """
    if chrom.startswith("chr"):
        chrom_tag = chrom[3:]  # Remove "chr" prefix
    else:
        chrom_tag = chrom

    if chrom_tag.lower() in ["m", "mitochondrion"]:
        return "ChrM"
    elif chrom_tag.upper() == "X":
        return "ChrX"
    elif chrom_tag.upper() == "Y":
        return "ChrY"
    else:
        return f"Chr{chrom_tag}"


def determine_allelic_state(chrom, gt):
    """
    Determines the allelic state based on the chromosome type and genotype.
    """
    if chrom == "ChrM":
        if gt in ["0/0", "0"]:
            return "homoplasmic"
    elif gt == "0/0":
        return "homozygous"
    elif gt == "0":
        return "hemizygous"
    return "unknown"


def parse_gvcf(gvcf_filename, build_version):
    """
    Parses the GVCF file and identifies contiguous studied regions.
    """
    studied_regions = []  # Contiguous regions (BED format)
    detailed_blocks = []  # Detailed blocks with FASTA reference allele info

    vcf_reader = vcf.Reader(filename=gvcf_filename)

    current_chromosome = None
    current_start = None
    current_end = None

    for record in vcf_reader:
        chrom = normalize_chromosome(record.CHROM)
        pos = record.POS
        end = record.INFO.get('END', pos)

        if current_start is None:
            current_chromosome = chrom
            current_start = pos
            current_end = end
        elif chrom == current_chromosome and pos <= current_end + 1:
            current_end = max(current_end, end)
        else:
            studied_regions.append({
                "CHROM": current_chromosome,
                "START": current_start - 1,
                "END": current_end
            })
            current_chromosome = chrom
            current_start = pos
            current_end = end

        gt = "0" if chrom == "ChrM" else record.samples[0].data.GT
        allelic_state = determine_allelic_state(chrom, gt)

        # Fetch reference allele from FASTA
        try:
            ref_allele = extract_chromosome_sequence(chrom.lstrip("Chr"), build_version)[pos - 1:end]
        except ValueError as e:
            print(f"Error extracting reference sequence for {chrom}:{pos}-{end}: {e}")
            ref_allele = ""

        detailed_blocks.append({
            "CHROM": chrom,
            "POS": pos,
            "END": end,
            "REF_ALLELE": ref_allele,
            "FILTER": record.FILTER,
            "GT": gt,
            "allelicState": allelic_state
        })

    if current_start is not None:
        studied_regions.append({
            "CHROM": current_chromosome,
            "START": current_start - 1,
            "END": current_end
        })

    return studied_regions, detailed_blocks


def generate_bed_json(studied_regions, build_version):
    """
    Generates BED JSON for contiguous studied regions.
    """
    bed_entries = []

    for region in studied_regions:
        chrom = region["CHROM"]
        start = region["START"]
        end = region["END"]

        try:
            chromosome_seq = extract_chromosome_sequence(chrom.lstrip("Chr"), build_version)
            ref_allele = chromosome_seq[start:end]
        except ValueError as e:
            print(f"Error processing region {chrom}:{start}-{end}: {e}")
            ref_allele = ""

        bed_entries.append({
            "CHROM": chrom,
            "START": start,
            "END": end,
            "REF_ALLELE": ref_allele
        })

    return bed_entries


def main(gvcf_file_path, build_version):
    """
    Main function to parse the GVCF file and generate BED JSON.
    """
    print(f"Parsing GVCF file: {gvcf_file_path} with build version {build_version}")
    studied_regions, detailed_blocks = parse_gvcf(gvcf_file_path, build_version)

    bed_json = generate_bed_json(studied_regions, build_version)

    with open("bed_output_with_ref_and_metadata.json", "w") as bed_file:
        json.dump(bed_json, bed_file, indent=2)
    print("BED JSON saved.")

    with open("non_variant_blocks.json", "w") as blocks_file:
        json.dump(detailed_blocks, blocks_file, indent=2)
    print("Detailed blocks saved.")


# Example usage
gvcf_file_path = "/Users/ethankakavetsis/Desktop/genomics-operations/utilities/NA18870.chr20.GRCh37.g.vcf"
build_version = 37

if __name__ == "__main__":
    main(gvcf_file_path, build_version)
