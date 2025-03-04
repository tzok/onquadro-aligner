#!/usr/bin/env python3
"""
DNA/RNA Sequence Alignment Tool

This program aligns two DNA or RNA sequences provided as command line arguments.
"""

import sys
import argparse


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Align two DNA/RNA sequences.")
    parser.add_argument("seq1", help="First DNA/RNA sequence")
    parser.add_argument("seq2", help="Second DNA/RNA sequence")
    return parser.parse_args()


def validate_sequence(sequence):
    """Validate if the input is a valid DNA/RNA sequence."""
    valid_chars = set("ATGCUN")
    if not all(char.upper() in valid_chars for char in sequence):
        return False
    return True


def calculate_score(seq1, seq2, i, j, consecutive_g_count):
    """
    Calculate alignment score with bonus for consecutive Gs.

    Args:
        seq1, seq2: The sequences being aligned
        i, j: Current positions in the sequences
        consecutive_g_count: Number of consecutive G matches so far

    Returns:
        Score for this position
    """
    # Base match/mismatch score
    if i < 0 or j < 0:
        return -1  # Gap penalty

    if seq1[i].upper() == seq2[j].upper():
        match_score = 2  # Basic match score

        # Bonus for G matches, especially consecutive ones
        if seq1[i].upper() == "G":
            match_score += 1  # Basic bonus for G
            match_score += consecutive_g_count  # Additional bonus for consecutive Gs

        return match_score
    else:
        return -1  # Mismatch penalty


def align_sequences(seq1, seq2):
    """
    Align two DNA/RNA sequences with emphasis on consecutive G matches.

    Uses dynamic programming to find the optimal alignment that
    maximizes the score, with special consideration for G nucleotides.
    """
    # Convert to uppercase for consistency
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Initialize the scoring matrix and traceback matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    traceback = [["" for _ in range(n + 1)] for _ in range(m + 1)]

    # Initialize consecutive G tracking matrix
    g_count_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # Initialize first row and column with gap penalties
    for i in range(1, m + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] - 1
        traceback[i][0] = "up"

    for j in range(1, n + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] - 1
        traceback[0][j] = "left"

    # Fill the matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Check if we have consecutive Gs
            consecutive_g_count = 0
            if (
                i > 1
                and j > 1
                and seq1[i - 2] == "G"
                and seq2[j - 2] == "G"
                and seq1[i - 1] == "G"
                and seq2[j - 1] == "G"
            ):
                consecutive_g_count = g_count_matrix[i - 1][j - 1]

            # Calculate scores for each possible move
            if seq1[i - 1] == "G" and seq2[j - 1] == "G":
                diagonal_score = score_matrix[i - 1][j - 1] + calculate_score(
                    seq1, seq2, i - 1, j - 1, consecutive_g_count
                )
                # Update G count for this position
                g_count_matrix[i][j] = g_count_matrix[i - 1][j - 1] + 1
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + calculate_score(
                    seq1, seq2, i - 1, j - 1, 0
                )
                g_count_matrix[i][j] = 0

            up_score = score_matrix[i - 1][j] - 1  # Gap in seq2
            left_score = score_matrix[i][j - 1] - 1  # Gap in seq1

            # Choose the best score
            if diagonal_score >= up_score and diagonal_score >= left_score:
                score_matrix[i][j] = diagonal_score
                traceback[i][j] = "diagonal"
            elif up_score >= left_score:
                score_matrix[i][j] = up_score
                traceback[i][j] = "up"
            else:
                score_matrix[i][j] = left_score
                traceback[i][j] = "left"

    # Traceback to find the alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = m, n

    while i > 0 or j > 0:
        if traceback[i][j] == "diagonal":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:  # 'left'
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse the sequences since we built them backwards
    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2


def display_alignment(aligned_seq1, aligned_seq2):
    """Display the aligned sequences with match indicators."""
    print("\nAlignment Result:")
    print(f"Sequence 1: {aligned_seq1}")

    # Create a match line to show matches between sequences
    match_line = ""
    for i in range(len(aligned_seq1)):
        if (
            i < len(aligned_seq2)
            and aligned_seq1[i] == aligned_seq2[i]
            and aligned_seq1[i] != "-"
        ):
            if aligned_seq1[i] == "G":
                match_line += "*"  # Special indicator for G matches
            else:
                match_line += "|"  # Regular match
        else:
            match_line += " "  # Mismatch or gap

    print(f"            {match_line}")
    print(f"Sequence 2: {aligned_seq2}")

    # Count and display consecutive G matches
    g_matches = 0
    consecutive_g_runs = []
    current_run = 0

    for i in range(min(len(aligned_seq1), len(aligned_seq2))):
        if aligned_seq1[i] == "G" and aligned_seq2[i] == "G":
            g_matches += 1
            current_run += 1
        elif current_run > 0:
            consecutive_g_runs.append(current_run)
            current_run = 0

    if current_run > 0:
        consecutive_g_runs.append(current_run)

    print(f"\nTotal G matches: {g_matches}")
    if consecutive_g_runs:
        print(f"Consecutive G runs: {consecutive_g_runs}")
        print(f"Longest consecutive G run: {max(consecutive_g_runs)}")


def main():
    """Main function to run the alignment tool."""
    args = parse_arguments()

    # Get sequences from command line arguments
    seq1 = args.seq1
    seq2 = args.seq2

    # Validate sequences
    if not validate_sequence(seq1) or not validate_sequence(seq2):
        print(
            "Error: Invalid sequence. Please use only A, T, G, C, U, or N characters."
        )
        sys.exit(1)

    # Print input sequences
    print(f"Input Sequence 1: {seq1}")
    print(f"Input Sequence 2: {seq2}")

    # Align sequences
    aligned_seq1, aligned_seq2 = align_sequences(seq1, seq2)

    # Display alignment
    display_alignment(aligned_seq1, aligned_seq2)


if __name__ == "__main__":
    main()
