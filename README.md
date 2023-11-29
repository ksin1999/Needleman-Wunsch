# NW
Dynamic Programming Project
# Needleman-Wunsch Sequence Alignment

This Python script was created as part of a Bioinformatics group project whereby the Needleman-Wunsch algorithm was implemented to solve a multiple sequence alignment problem. 
This script takes two input sequences and calculates the optimal alignment and alignment score. To 

## Features

- **Global Sequence Alignment:** Uses the Needleman-Wunsch algorithm to align two sequences globally.
- **Scoring System:** Allows customization of match, mismatch, and gap scores.
- **Traceback Visualization:** Provides a visual representation of the alignment.

# Example Usage
# Dictionary of sequence pairs
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

Customization
You can customize the scoring system by modifying the match, mismatch, and gap parameters in the script.

python
Copy code
# Customization
match = 1
mismatch = -1
gap = -2
