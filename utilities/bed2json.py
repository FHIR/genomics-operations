import pyranges
from collections import OrderedDict


def bed2json(bed_filename=None):
    output_json = OrderedDict()
    output_json["BedID"] = bed_filename[5:]
    if not (bed_filename):
        raise Exception('You must provide bed_filename')

    cr = []
    df = pyranges.read_bed(bed_filename, as_df=True)
    df = df.to_dict('list')
    for i in range(len(df['Chromosome'])):
        temp = {}
        temp['Chromosome'] = df['Chromosome'][i]
        temp['Start'] = df['Start'][i]
        temp['End'] = df['End'][i]
        cr.append(temp)

    output_json["BED"] = cr
    return output_json
