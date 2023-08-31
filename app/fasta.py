from pathlib import Path
from urllib.request import urlretrieve

FASTA_INDEX_RELEASE_TAG = '113c119'


def download_fasta():
    try:
        # Make sure the parent folder exists
        Path('FASTA').mkdir(exist_ok=True)

        for build in ['GRCh37', 'GRCh38']:
            filename = build + '_latest_genomic.fna.gz'
            filepath = 'FASTA/' + filename

            # Download FASTA files
            if not Path(filepath).is_file():
                urlretrieve('https://github.com/FHIR/genomics-operations/releases/download/' + FASTA_INDEX_RELEASE_TAG + '/' + filename, filepath)

            indexFilename = filename + '.fxi'
            indexFilepath = 'FASTA/' + indexFilename

            # Download FASTA index files
            if not Path(indexFilepath).is_file():
                urlretrieve('https://github.com/FHIR/genomics-operations/releases/download/' + FASTA_INDEX_RELEASE_TAG + '/' + indexFilename, indexFilepath)

    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
