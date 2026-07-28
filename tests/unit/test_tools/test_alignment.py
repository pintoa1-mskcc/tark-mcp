from tark_mcp.tools.alignment import needleman_wunsch


def test_identical_sequences_align_with_no_gaps():
    aligned_a, aligned_b = needleman_wunsch("MATS", "MATS")
    assert aligned_a == "MATS"
    assert aligned_b == "MATS"


def test_single_substitution_aligns_without_gaps():
    aligned_a, aligned_b = needleman_wunsch("MAT", "MST")
    assert aligned_a == "MAT"
    assert aligned_b == "MST"


def test_insertion_in_second_sequence_produces_gap_in_first():
    aligned_a, aligned_b = needleman_wunsch("AC", "ABC")
    assert aligned_a == "A-C"
    assert aligned_b == "ABC"


def test_deletion_in_second_sequence_produces_gap_in_second():
    aligned_a, aligned_b = needleman_wunsch("ABC", "AC")
    assert aligned_a == "ABC"
    assert aligned_b == "A-C"


def test_aligned_sequences_are_always_equal_length():
    aligned_a, aligned_b = needleman_wunsch("MPIGSKERP", "MVLSPAD")
    assert len(aligned_a) == len(aligned_b)


def test_empty_sequence_aligns_as_all_gaps():
    aligned_a, aligned_b = needleman_wunsch("", "ABC")
    assert aligned_a == "---"
    assert aligned_b == "ABC"
