import os
from os.path import isdir

import pyard

# Make sure the pyard folder exists locally
if not isdir('./data/pyard'):
    exit("Missing pyard folder. Please run fetch_utilities_data.sh!")

pyard_database_version = os.getenv('PYARD_DATABASE_VERSION', '3580')
ard = pyard.init(data_dir='./data/pyard', cache_size=1, imgt_version=pyard_database_version)


def redux(glstring, redux_type="lgx"):
    return ard.redux(glstring, redux_type)
