import pyfastx
from threading import Lock

# Fasta file handles cache
fasta_cache = {}
fasta_lock = Lock()

# File paths
BUILD37_FILE = 'FASTA/GCF_000001405.25_GRCh37.p13_genomic.fna'
BUILD38_FILE = 'FASTA/GCF_000001405.40_GRCh38.p14_genomic.fna'

# Dictionaries for GRCh37 and GRCh38 chromosomes
grch37_chromosomes = {
    "1": "NC_000001.10", "2": "NC_000002.11", "3": "NC_000003.11", "4": "NC_000004.11", "5": "NC_000005.9",
    "6": "NC_000006.11", "7": "NC_000007.13", "8": "NC_000008.10", "9": "NC_000009.11", "10": "NC_000010.10",
    "11": "NC_000011.9", "12": "NC_000012.11", "13": "NC_000013.10", "14": "NC_000014.8", "15": "NC_000015.9",
    "16": "NC_000016.9", "17": "NC_000017.10", "18": "NC_000018.9", "19": "NC_000019.9", "20": "NC_000020.10",
    "21": "NC_000021.8", "22": "NC_000022.10", "X": "NC_000023.10", "Y": "NC_000024.9", "MIT": "NC_012920.1"
}
grch38_chromosomes = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12", "5": "NC_000005.10",
    "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12", "10": "NC_000010.11",
    "11": "NC_000011.10", "12": "NC_000012.12", "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10", "MIT": "NC_012920.1"
}

def get_fasta(file):
    """
    Retrieves the pyfastx.Fasta object for a given file, using a caching mechanism.

    Args:
        file (str): Path to the FASTA file.

    Returns:
        pyfastx.Fasta: An indexed FASTA object.
    """
    with fasta_lock:
        if file not in fasta_cache:
            try:
                fasta = pyfastx.Fasta(file)
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                raise
            fasta_cache[file] = fasta
        return fasta_cache[file]

def extract_chromosome_sequence(chromosome_tag, file_version):
    """
    Extracts the sequence of the specified chromosome using a cached FASTA file.

    Args:
        chromosome_tag (str): Chromosome tag (1-22, X, Y, or MIT).
        file_version (int): Version of the genome build (37 for GRCh37, 38 for GRCh38).

    Returns:
        str: The full sequence of the specified chromosome.
    """
    if file_version == 37:
        file_path = BUILD37_FILE
        chromosome = grch37_chromosomes.get(chromosome_tag)
    elif file_version == 38:
        file_path = BUILD38_FILE
        chromosome = grch38_chromosomes.get(chromosome_tag)
    else:
        raise ValueError("Invalid file version. Use 37 for GRCh37 or 38 for GRCh38.")

    if not chromosome:
        raise ValueError(f"Invalid chromosome tag: {chromosome_tag}. Must be 1-22, X, Y, or MIT.")

    fasta = get_fasta(file_path)
    if chromosome in fasta:
        return fasta[chromosome].seq
    else:
        raise ValueError(f"Chromosome {chromosome} not found in {file_path}")

def main():
    # Example: Extract chromosome 20 and mitochondrion sequences
    chromosome_tag = "20"  # Chromosome 20
    mito_tag = "MIT"       # Mitochondrion

    sequence_37_chr20 = extract_chromosome_sequence(chromosome_tag, 37)
    sequence_38_chr20 = extract_chromosome_sequence(chromosome_tag, 38)

    sequence_37_mito = extract_chromosome_sequence(mito_tag, 37)
    sequence_38_mito = extract_chromosome_sequence(mito_tag, 38)

    # Print lengths and previews
    print(f"\nLength of Chromosome {chromosome_tag} (GRCh37): {len(sequence_37_chr20)}")
    print(f"Preview (GRCh37): {sequence_37_chr20[:100]}...")

    print(f"\nLength of Chromosome {chromosome_tag} (GRCh38): {len(sequence_38_chr20)}")
    print(f"Preview (GRCh38): {sequence_38_chr20[:100]}...")

    print(f"\nLength of Mitochondrion (GRCh37): {len(sequence_37_mito)}")
    print(f"Preview (GRCh37): {sequence_37_mito[:100]}...")

    print(f"\nLength of Mitochondrion (GRCh38): {len(sequence_38_mito)}")
    print(f"Preview (GRCh38): {sequence_38_mito[:100]}...")

if __name__ == "__main__":
    main()
