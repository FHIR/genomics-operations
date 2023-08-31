from pathlib import Path
from urllib.request import urlretrieve
from utilities.SPDI_Normalization import init_fasta


def download_fasta():
    try:
        # Make sure the parent folder exists
        Path('FASTA').mkdir(exist_ok=True)

        for build in ['GRCh37', 'GRCh38']:
            filename = build + '_latest_genomic.fna.gz'
            filepath = 'FASTA/' + filename

            # Download files
            if not Path(filepath).is_file():
                urlretrieve('https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/' + build + '_latest/refseq_identifiers/' + filename, filepath)

        # Build indexes
        init_fasta()
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
