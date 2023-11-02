import pyfastx
from threading import Lock


# Fasta file handles cache
fasta_cache = {}
fasta_lock = Lock()

BUILD37_FILE = 'FASTA/GRCh37_latest_genomic.fna.gz'
BUILD38_FILE = 'FASTA/GRCh38_latest_genomic.fna.gz'


def get_fasta(file):
    with fasta_lock:
        if file not in fasta_cache:
            try:
                fasta = pyfastx.Fasta(file)
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                raise
            fasta_cache[file] = fasta
        return fasta_cache[file]


def get_normalized_spdi(ref_seq, pos, ref, alt, build):
    if build == 'GRCh37':
        GRCh37_ref_seq_fasta = get_fasta(BUILD37_FILE)
        ref_seq_fasta = GRCh37_ref_seq_fasta[ref_seq]
    elif build == 'GRCh38':
        GRCh38_ref_seq_fasta = get_fasta(BUILD38_FILE)
        ref_seq_fasta = GRCh38_ref_seq_fasta[ref_seq]

    # Step 0
    start = pos
    end = pos + len(ref)

    # Step 1
    # a. Trim Common Suffix

    small = len(alt) if len(alt) <= len(ref) else len(ref)

    suffix_length = 0

    last = -1
    while small != 0:
        if ref[last] == alt[last]:
            suffix_length += 1
            end -= 1
            last -= 1
            small -= 1
        else:
            break

    if suffix_length > 0:
        ref = ref[0:-suffix_length]
        alt = alt[0:-suffix_length]

        # b. Trim Common Prefix

    small = len(alt) if len(alt) <= len(ref) else len(ref)

    prefix_length = 0

    first = 0
    while small != 0:
        if ref[first] == alt[first]:
            prefix_length += 1
            start += 1
            first += 1
            small -= 1
        else:
            break

    if prefix_length > 0:
        ref = ref[prefix_length:]
        alt = alt[prefix_length:]

    # Step 2
        # a. Check if both are empty
        # Not Possible in this case

        # b. Check if both are non-empty

    if ref and alt:
        return f"{ref_seq}:{start}:{ref}:{alt}"

        # c. Check if one is empty

    if alt:
        allele = alt
        type_ = 'I'
    else:
        allele = ref
        type_ = 'D'

    # Step 3
        # a. Left Roll

    if type_ == 'I':
        left_roll_bound = start

        temp_allele = allele

        while (str(ref_seq_fasta[left_roll_bound-1]).upper() == temp_allele[-1].upper()):
            temp_allele = str(ref_seq_fasta[left_roll_bound-1]).upper() + temp_allele[:-1]
            left_roll_bound -= 1

            # . b. Right Roll

        right_roll_bound = start

        temp_allele = allele

        while (str(ref_seq_fasta[right_roll_bound]).upper() == temp_allele[0].upper()):
            temp_allele = temp_allele[1:] + str(ref_seq_fasta[right_roll_bound]).upper()
            right_roll_bound += 1

        # Step 4
            # a. Prepend
        if ref_seq_fasta[left_roll_bound:start]:
            ref = str(ref_seq_fasta[left_roll_bound:start]) + ref
            alt = str(ref_seq_fasta[left_roll_bound:start]) + alt

            # b. append
        if ref_seq_fasta[start:right_roll_bound]:
            ref = ref + str(ref_seq_fasta[start:right_roll_bound])
            alt = alt + str(ref_seq_fasta[start:right_roll_bound])

            # c. Modify Start and End, and Return SPDI
        start = left_roll_bound
        end = right_roll_bound

        return (f"{ref_seq}:{start}:{ref.upper()}:{alt.upper()}")

    else:
        left_roll_bound = start

        while (str(ref_seq_fasta[left_roll_bound-1]).upper() == str(ref_seq_fasta[left_roll_bound + (len(allele)-1)]).upper()):
            left_roll_bound -= 1

            # . b. Right Roll

        right_roll_bound = start

        while (str(ref_seq_fasta[right_roll_bound]).upper() == (ref_seq_fasta[right_roll_bound + (len(allele))]).upper()):
            right_roll_bound += 1

        # Step 4
            # a. Prepend
        if ref_seq_fasta[left_roll_bound:start]:
            ref = str(ref_seq_fasta[left_roll_bound:start]) + ref
            alt = str(ref_seq_fasta[left_roll_bound:start]) + alt

            # b. append
        if ref_seq_fasta[start:right_roll_bound]:
            ref = ref + str(ref_seq_fasta[start:right_roll_bound])
            alt = alt + str(ref_seq_fasta[start:right_roll_bound])

            # c. Modify Start and End, and Return SPDI
        start = left_roll_bound
        end = right_roll_bound

        return (f"{ref_seq}:{start}:{ref.upper()}:{alt.upper()}")
