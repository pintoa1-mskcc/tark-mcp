from __future__ import annotations


def needleman_wunsch(
    seq_a: str,
    seq_b: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> tuple[str, str]:
    """Global (Needleman-Wunsch) alignment of two sequences.

    Returns a pair of equal-length strings representing the alignment of
    seq_a and seq_b, with '-' inserted at gap positions. Uses a linear gap
    penalty and a simple match/mismatch scoring scheme (no substitution
    matrix).
    """
    n, m = len(seq_a), len(seq_b)

    # score[i][j] = best alignment score for seq_a[:i] vs seq_b[:j]
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * gap
    for j in range(1, m + 1):
        score[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            score[i][j] = max(
                score[i - 1][j - 1] + s,
                score[i - 1][j] + gap,
                score[i][j - 1] + gap,
            )

    aligned_a: list[str] = []
    aligned_b: list[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            if score[i][j] == score[i - 1][j - 1] + s:
                aligned_a.append(seq_a[i - 1])
                aligned_b.append(seq_b[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and score[i][j] == score[i - 1][j] + gap:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
            continue
        aligned_a.append("-")
        aligned_b.append(seq_b[j - 1])
        j -= 1

    aligned_a.reverse()
    aligned_b.reverse()
    return "".join(aligned_a), "".join(aligned_b)
