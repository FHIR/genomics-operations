from pathlib import Path
from pyfastx import Fasta
from urllib.request import urlretrieve


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
            if not Path(filepath + '.fxi').is_file():
                Fasta(filepath)
    except Exception as error:
        print(error)
