import json
import vcf
from tqdm import tqdm
from common import extract_chromosome_sequence


def parse_gvcf(gvcf_filename, build_version):
    """
    Parses the GVCF file and identifies contiguous studied regions.
    Ensures REF_ALLELE is converted to uppercase and handles multi-base REF alleles correctly.
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
        end = record.INFO.get('END', pos)  # Get the END info, or use POS if not available
        ref_len_end = len(record.REF) + pos  # Length of the REF allele + POS to calculate the true end

        # Check if the current record is part of the same contiguous region
        if current_start is None:
            current_chromosome = chrom
            current_start = pos
            current_end = max(end, ref_len_end)  # Ensure the region ends at the correct position
        elif chrom == current_chromosome and pos <= current_end + 1:
            # The variants are contiguous, so extend the region
            current_end = max(current_end, end, ref_len_end)
        else:
            # Save the previous region before starting a new one
            studied_regions.append({
                "CHROM": current_chromosome,
                "START": current_start - 1,
                "END": current_end
            })
            # Start a new region
            current_chromosome = chrom
            current_start = pos
            current_end = max(end, ref_len_end)

        gt = "0" if chrom == "chrM" else record.samples[0].data.GT
        allelic_state = determine_allelic_state(chrom, gt)

        # Fetch reference allele from FASTA and capitalize
        try:
            ref_allele = extract_chromosome_sequence(chrom.lstrip("chr"), build_version)[pos - 1:end].upper()
        except ValueError as e:
            print(f"Error extracting reference sequence for {chrom}:{pos}-{end}: {e}")
            ref_allele = ""

        detailed_blocks.append({
            "CHROM": chrom,
            "POS": pos,
            "END": end,
            "REF_ALLELE": ref_allele,  # Capitalized REF_ALLELE
            "FILTER": record.FILTER,
            "GT": gt,
            "allelicState": allelic_state
        })

    # Append the last region
    if current_start is not None:
        studied_regions.append({
            "CHROM": current_chromosome,
            "START": current_start - 1,
            "END": current_end
        })

    return studied_regions, detailed_blocks



def normalize_chromosome(chrom):
    """
    Normalizes chromosome names to the `Chr` format.
    """
    if chrom.startswith("chr"):
        chrom_tag = chrom[3:]  # Remove "chr" prefix
    else:
        chrom_tag = chrom

    if chrom_tag.lower() in ["m", "mitochondrion"]:
        return "chrM"
    elif chrom_tag.upper() == "X":
        return "chrX"
    elif chrom_tag.upper() == "Y":
        return "chrY"
    else:
        return f"chr{chrom_tag}"


def determine_allelic_state(chrom, gt):
    """
    Determines the allelic state based on the chromosome type and genotype.
    """
    if chrom == "chrM":
        if gt in ["0/0", "0"]:
            return "homoplasmic"
    elif gt == "0/0":
        return "homozygous"
    elif gt == "0":
        return "hemizygous"
    return "unknown"



