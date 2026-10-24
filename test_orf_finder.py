from orf_finder import find_orfs


def test_simple_forward_orf():
    orfs = find_orfs("ATGAAATAG", min_aa=1)
    assert len(orfs) == 1
    assert orfs[0]["protein"] == "MK"
    assert orfs[0]["strand"] == "+"
    assert (orfs[0]["start"], orfs[0]["end"]) == (1, 9)


def test_reverse_strand_orf():
    # reverse complement of ATGAAATAG
    orfs = find_orfs("CTATTTCAT", min_aa=1)
    assert len(orfs) == 1
    assert orfs[0]["strand"] == "-"
    assert orfs[0]["protein"] == "MK"
    assert (orfs[0]["start"], orfs[0]["end"]) == (1, 9)


def test_min_length_filter():
    assert find_orfs("ATGAAATAG", min_aa=10) == []


def test_no_stop_codon_is_ignored():
    assert find_orfs("ATGAAAAAA", min_aa=1) == []