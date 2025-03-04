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
    parser.add_argument("-n", "--num-alignments", type=int, default=5,
                        help="Number of alignments to generate (default: 5)")
    parser.add_argument("-s", "--score-threshold", type=float, default=0.8,
                        help="Score threshold as a fraction of the optimal score (0.0-1.0, default: 0.8)")
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


def align_sequences(seq1, seq2, num_alignments=5, score_threshold=0.8):
    """
    Align two DNA/RNA sequences with emphasis on consecutive G matches.
    
    Uses dynamic programming to find multiple optimal and suboptimal alignments
    that maximize the score, with special consideration for G nucleotides.
    
    Args:
        seq1, seq2: The sequences to align
        num_alignments: Maximum number of alignments to return
        score_threshold: Minimum score threshold as a fraction of optimal score (0.0-1.0)
        
    Returns:
        List of tuples (aligned_seq1, aligned_seq2, score) sorted by score
    """
    # Convert to uppercase for consistency
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Initialize the scoring matrix and traceback matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Store all possible moves at each cell instead of just the best one
    traceback = [[[] for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Initialize consecutive G tracking matrix
    g_count_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # Initialize first row and column with gap penalties
    for i in range(1, m + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] - 1
        traceback[i][0].append(("up", score_matrix[i][0]))

    for j in range(1, n + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] - 1
        traceback[0][j].append(("left", score_matrix[0][j]))

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

            # Store all possible moves and their scores
            moves = []
            
            # Find the best score
            best_score = max(diagonal_score, up_score, left_score)
            score_matrix[i][j] = best_score
            
            # Add all moves that are within threshold of the best score
            if diagonal_score >= best_score * score_threshold:
                moves.append(("diagonal", diagonal_score))
            if up_score >= best_score * score_threshold:
                moves.append(("up", up_score))
            if left_score >= best_score * score_threshold:
                moves.append(("left", left_score))
                
            # Sort moves by score (highest first)
            moves.sort(key=lambda x: x[1], reverse=True)
            traceback[i][j] = moves

    # Generate multiple alignments using backtracking
    alignments = []
    optimal_score = score_matrix[m][n]
    min_score_threshold = optimal_score * score_threshold
    
    # Use depth-first search to find multiple paths
    def backtrack(i, j, aligned1, aligned2, current_score):
        if i == 0 and j == 0:
            # We've reached the beginning of both sequences
            alignments.append((
                "".join(reversed(aligned1)), 
                "".join(reversed(aligned2)), 
                current_score
            ))
            return
            
        # Try all possible moves from this cell
        for move, move_score in traceback[i][j]:
            if move == "diagonal":
                backtrack(
                    i - 1, 
                    j - 1, 
                    aligned1 + [seq1[i - 1]], 
                    aligned2 + [seq2[j - 1]], 
                    current_score
                )
            elif move == "up":
                backtrack(
                    i - 1, 
                    j, 
                    aligned1 + [seq1[i - 1]], 
                    aligned2 + ["-"], 
                    current_score
                )
            elif move == "left":
                backtrack(
                    i, 
                    j - 1, 
                    aligned1 + ["-"], 
                    aligned2 + [seq2[j - 1]], 
                    current_score
                )
                
            # Stop if we have enough alignments
            if len(alignments) >= num_alignments:
                return
    
    # Start backtracking from the bottom-right cell
    backtrack(m, n, [], [], optimal_score)
    
    # Sort alignments by score (highest first) and remove duplicates
    unique_alignments = []
    seen = set()
    
    for aligned_seq1, aligned_seq2, score in sorted(alignments, key=lambda x: x[2], reverse=True):
        alignment_key = (aligned_seq1, aligned_seq2)
        if alignment_key not in seen and score >= min_score_threshold:
            seen.add(alignment_key)
            unique_alignments.append((aligned_seq1, aligned_seq2, score))
            
            # Stop if we have enough unique alignments
            if len(unique_alignments) >= num_alignments:
                break
    
    return unique_alignments


def display_alignment(aligned_seq1, aligned_seq2, score=None, alignment_num=None):
    """
    Display a single alignment with match indicators and statistics.
    
    Args:
        aligned_seq1, aligned_seq2: The aligned sequences
        score: The alignment score (optional)
        alignment_num: The alignment number (optional)
    """
    header = "\nAlignment"
    if alignment_num is not None:
        header += f" #{alignment_num}"
    if score is not None:
        header += f" (Score: {score})"
    print(header + ":")
    
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

    print(f"G matches: {g_matches}")
    if consecutive_g_runs:
        print(f"Consecutive G runs: {consecutive_g_runs}")
        print(f"Longest consecutive G run: {max(consecutive_g_runs) if consecutive_g_runs else 0}")


def display_all_alignments(alignments):
    """
    Display multiple alignments with their scores.
    
    Args:
        alignments: List of tuples (aligned_seq1, aligned_seq2, score)
    """
    print(f"\nFound {len(alignments)} alignment(s):")
    
    for i, (aligned_seq1, aligned_seq2, score) in enumerate(alignments, 1):
        display_alignment(aligned_seq1, aligned_seq2, score, i)


def main():
    """Main function to run the alignment tool."""
    args = parse_arguments()

    # Get sequences from command line arguments
    seq1 = args.seq1
    seq2 = args.seq2
    num_alignments = args.num_alignments
    score_threshold = args.score_threshold

    # Validate sequences
    if not validate_sequence(seq1) or not validate_sequence(seq2):
        print(
            "Error: Invalid sequence. Please use only A, T, G, C, U, or N characters."
        )
        sys.exit(1)

    # Print input sequences
    print(f"Input Sequence 1: {seq1}")
    print(f"Input Sequence 2: {seq2}")
    print(f"Generating up to {num_alignments} alignments with score threshold {score_threshold}")

    # Align sequences and get multiple alignments
    alignments = align_sequences(seq1, seq2, num_alignments, score_threshold)

    # Display all alignments
    display_all_alignments(alignments)


if __name__ == "__main__":
    main()
