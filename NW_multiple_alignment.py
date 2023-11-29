import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Perform global sequence alignment using the Needleman-Wunsch algorithm.

    Parameters:
    - seq1 (str): The first input sequence.
    - seq2 (str): The second input sequence.
    - match (int): Score for a match.
    - mismatch (int): Score for a mismatch.
    - gap (int): Score for a gap.

    Returns:
    - Tuple: Aligned sequences and the alignment score.
    """
    # Create matrices
    score_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1))
    match_matrix = np.zeros((len(seq1), len(seq2)))

    # Fill match/mismatch matrix
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            match_matrix[i][j] = match if seq1[i] == seq2[j] else mismatch

    # Initialization
    for i in range(len(seq1) + 1):
        score_matrix[i][0] = i * gap
    for j in range(len(seq2) + 1):
        score_matrix[0][j] = j * gap

    # Matrix filling
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score_matrix[i][j] = max(
                score_matrix[i - 1][j - 1] + match_matrix[i - 1][j - 1],
                score_matrix[i - 1][j] + gap,
                score_matrix[i][j - 1] + gap
            )

    # Traceback
    aligned_seq1, aligned_seq2 = "", ""
    ti, tj = len(seq1), len(seq2)

    while ti > 0 or tj > 0:
        diag = score_matrix[ti - 1][tj - 1]
        left = score_matrix[ti][tj - 1]
        top = score_matrix[ti - 1][tj]

        if diag != left + gap != top + gap:
            if diag == max(diag, left + gap, top + gap):
                aligned_seq1, aligned_seq2, tj, ti = seq1[ti - 1] + aligned_seq1, seq2[tj - 1] + aligned_seq2, tj - 1, ti - 1
            elif left + gap == max(diag, left + gap, top + gap):
                aligned_seq1, aligned_seq2, tj = "-" + aligned_seq1, seq2[tj - 1] + aligned_seq2, tj - 1
            else:
                aligned_seq1, aligned_seq2, ti = seq1[ti - 1] + aligned_seq1, "-" + aligned_seq2, ti - 1
        else:
            opt = max(diag, left, top)
            if opt == left:
                aligned_seq1, aligned_seq2, tj = "-" + aligned_seq1, seq2[tj - 1] + aligned_seq2, tj - 1
            elif opt == diag:
                aligned_seq1, aligned_seq2, tj, ti = seq1[ti - 1] + aligned_seq1, seq2[tj - 1] + aligned_seq2, tj - 1, ti - 1
            else:
                aligned_seq1, aligned_seq2, ti = seq1[ti - 1] + aligned_seq1, "-" + aligned_seq2, ti - 1

    # Return results
    alignment_score = score_matrix[len(seq1)][len(seq2)]
    return aligned_seq1, aligned_seq2, alignment_score

# Example usage
seq_dict = {
    "pseq_1": ["GATTACA", "ATTACATTAC"],
    "pseq_2": ["ATATATATATA", "ATATATATAT"],
    "pseq_3": ["AATAATAATAAT", "AAAATAAATAAA"],
    "pseq_4": ["ATATACACACA", "ATATGTATACAT"]
}

selected_seq = "pseq_1"
seq1, seq2 = seq_dict[selected_seq]

aligned_seq1, aligned_seq2, alignment_score = needleman_wunsch(seq1, seq2)

# Display results
print(f"Alignment 1: {aligned_seq1}")
print(f"Alignment 2: {aligned_seq2}")
print(f"Alignment Score: {alignment_score}")