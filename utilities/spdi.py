from bioutils.normalize import NormalizationMode
from bioutils.normalize import normalize as spdi_normalize


def normalize(ref_seq, acc, pos, ref, alt):
    new_ival, new_alleles = spdi_normalize(sequence=ref_seq,
                                           interval=(pos, pos + len(ref)),
                                           alleles=(None, alt),
                                           mode=NormalizationMode.EXPAND,
                                           anchor_length=0,
                                           trim=ref != alt)
    return (f"{acc}:{new_ival[0]}:{new_alleles[0]}:{new_alleles[1]}")
