#!/bin/sh

set -eu

. ./.env

if [ -d ./data ]; then
    echo "files already fetched."
    exit 0
fi

mkdir -p ./data/pyard
(
    cd ./data/pyard

    echo "Downloading py-ard database..."

    curl -sLO https://github.com/FHIR/genomics-operations/releases/download/${UTILITIES_DATA_VERSION}/pyard.sqlite3.tar.gz

    echo "Extracting py-ard database..."

    tar -xzf pyard.sqlite3.tar.gz

    echo "Finished extracting py-ard database."
)
