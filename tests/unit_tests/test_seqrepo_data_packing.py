import pytest

from utilities.pack_utilities_data import pack_ref_seq
from utilities.spdi_refseq import RefSeq


@pytest.mark.parametrize(
    "name, ref_seq, exp_packed_ref_seq",
    [
        ('happy flow', 'ACGT', b'\x88\xf6'),
        ('empty RefSeq', '', b''),
        ('2 nucleotides -> one byte', 'GT', b'\xda'),
        ('3 nucleotides -> two bytes', 'GTA', b'\x1a\xfe'),
        ('5 nucleotides -> two bytes', 'GTACT', b'\x1a\xb2'),
        ('6 nucleotides -> three bytes', 'GTACTG', b'\x1a\x32\xfd'),
        ('8 nucleotides -> three bytes', 'GTACTGAT', b'\x1a\x32\x61'),
        ('9 nucleotides -> four bytes', 'GTACTGATC', b'\x1a\x32\x61\xf9'),
        ('all other codes get mapped to N (including invalid ones)', 'FO0b4R', b'\xff\xff\xff')
    ],
)
def test_pack_ref_seq(name, ref_seq, exp_packed_ref_seq):
    assert pack_ref_seq(ref_seq) == exp_packed_ref_seq, name


@pytest.mark.parametrize(
    "name, packed_ref_seq, exp_ref_seq",
    [
        ('happy flow', b'\x88\xf6', 'ACGT'),
        ('empty RefSeq', b'', ''),
        ('2 nucleotides -> one byte', b'\xda', 'GT'),
        ('3 nucleotides -> two bytes', b'\x1a\xfe', 'GTA'),
        ('5 nucleotides -> two bytes', b'\x1a\xb2', 'GTACT'),
        ('6 nucleotides -> three bytes', b'\x1a\x32\xfd', 'GTACTG'),
        ('8 nucleotides -> three bytes', b'\x1a\x32\x61', 'GTACTGAT'),
        ('9 nucleotides -> four bytes', b'\x1a\x32\x61\xf9', 'GTACTGATC'),
        ('Ns are decoded correctly', b'\xff\xff\xff', 'NNNNNN')
    ],
)
def test_unpack_ref_seq(name, packed_ref_seq, exp_ref_seq):
    assert RefSeq(packed_ref_seq, len(exp_ref_seq))[:] == exp_ref_seq, name


def test_unpack_ref_seq_error():
    with pytest.raises(NotImplementedError) as exc_info:
        RefSeq(b'\x06', 1)[:]
    assert str(exc_info.value) == "unexpected hex-encoded code: '6'"
