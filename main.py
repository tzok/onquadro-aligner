#!/usr/bin/env python3
"""
DNA/RNA Sequence Alignment Tool

This program aligns two DNA or RNA sequences provided as command line arguments.
"""

import sys
import argparse
import json
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class QuadruplexDotBracket:
    """Class for storing quadruplex dot-bracket notation data."""

    sequence: str
    structure: str
    chi: str
    loop: str

    def __str__(self):
        """String representation of the quadruplex."""
        return (
            f"Sequence: {self.sequence}\n"
            f"Structure: {self.structure}\n"
            f"Chi: {self.chi}\n"
            f"Loop: {self.loop}"
        )

    def validate(self):
        """
        Validate that all fields have the same length.

        Returns:
            bool: True if valid, False otherwise
        """
        # Check if all fields have the same length
        fields = [self.sequence, self.structure, self.chi, self.loop]
        if len(set(len(field) for field in fields)) != 1:
            return False

        return True

    def get_segments(self):
        """
        Get the segments of the quadruplex (parts separated by '-').

        Returns:
            List of tuples (sequence_segment, structure_segment, chi_segment, loop_segment)
        """
        seq_segments = self.sequence.split("-")
        struct_segments = self.structure.split("-")
        chi_segments = self.chi.split("-")
        loop_segments = self.loop.split("-")

        return list(zip(seq_segments, struct_segments, chi_segments, loop_segments))


def read_quadruplex_json(file_path):
    """
    Read a JSON file containing quadruplexDotBracket objects.

    Args:
        file_path: Path to the JSON file

    Returns:
        List of QuadruplexDotBracket objects
    """
    try:
        with open(file_path, "r") as f:
            data = json.load(f)

        quadruplexes = []

        # Handle different JSON structures
        if isinstance(data, dict):
            # Check if the data contains a "quadruplexDotBracket" object
            if "quadruplexDotBracket" in data:
                quad_data = data["quadruplexDotBracket"]
                quad = parse_quadruplex_object(quad_data)
                if quad and quad.validate():
                    quadruplexes.append(quad)
                else:
                    print(
                        f"Warning: Invalid quadruplexDotBracket object found in {file_path} (fields have different lengths or substring structure)"
                    )
            else:
                # Try to parse as a direct quadruplex object
                quad = parse_quadruplex_object(data)
                if quad and quad.validate():
                    quadruplexes.append(quad)
                else:
                    print(
                        f"Warning: Invalid quadruplex object found in {file_path} (fields have different lengths or substring structure)"
                    )
        elif isinstance(data, list):
            # Array of objects - could be array of quadruplexDotBracket containers or direct objects
            for item in data:
                if isinstance(item, dict) and "quadruplexDotBracket" in item:
                    quad_data = item["quadruplexDotBracket"]
                    quad = parse_quadruplex_object(quad_data)
                else:
                    quad = parse_quadruplex_object(item)

                if quad and quad.validate():
                    quadruplexes.append(quad)
                else:
                    print(
                        f"Warning: Invalid quadruplex object found in {file_path} (fields have different lengths or substring structure)"
                    )

        return quadruplexes

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return []
    except json.JSONDecodeError:
        print(f"Error: File '{file_path}' contains invalid JSON.")
        return []
    except Exception as e:
        print(f"Error reading quadruplex data: {str(e)}")
        return []


def parse_quadruplex_object(data):
    """
    Parse a single quadruplexDotBracket object from JSON data.

    Args:
        data: Dictionary containing quadruplex data

    Returns:
        QuadruplexDotBracket object or None if parsing fails
    """
    try:
        # Extract required fields
        sequence = str(data.get("sequence", ""))
        structure = str(data.get("structure", ""))

        # Handle chi field - ensure it's a string
        chi_data = data.get("chi", "")
        if not isinstance(chi_data, str):
            chi = str(chi_data)
        else:
            chi = chi_data

        # Handle loop field - ensure it's a string
        loop_data = data.get("loop", "")
        if not isinstance(loop_data, str):
            if isinstance(loop_data, list):
                loop = "-".join(str(x) for x in loop_data)
            else:
                loop = str(loop_data)
        else:
            loop = loop_data

        # Create and return the object
        return QuadruplexDotBracket(sequence, structure, chi, loop)

    except Exception as e:
        print(f"Error parsing quadruplex object: {str(e)}")
        return None


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Align DNA/RNA sequences against quadruplex structures."
    )
    parser.add_argument("sequence", help="DNA/RNA sequence to align")
    parser.add_argument(
        "-d",
        "--directory",
        help="Directory containing JSON files with quadruplex data",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--num-alignments",
        type=int,
        default=5,
        help="Number of alignments to generate per quadruplex (default: 5)",
    )
    parser.add_argument(
        "-s",
        "--score-threshold",
        type=float,
        default=0.8,
        help="Score threshold as a fraction of the optimal score (0.0-1.0, default: 0.8)",
    )
    parser.add_argument(
        "-t",
        "--top-results",
        type=int,
        default=10,
        help="Number of top-scoring alignments to display (default: 10)",
    )
    return parser.parse_args()


def validate_sequence(sequence):
    """Validate if the input is a valid DNA/RNA sequence."""
    valid_chars = set("ATGCUN-")  # Include '-' as a valid character
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


def align_sequences(seq1, seq2, num_alignments=5, score_threshold=0.8, pre_computed_matrix=None):
    """
    Align two DNA/RNA sequences with emphasis on consecutive G matches.

    Uses dynamic programming to find multiple optimal and suboptimal alignments
    that maximize the score, with special consideration for G nucleotides.

    Args:
        seq1, seq2: The sequences to align
        num_alignments: Maximum number of alignments to return
        score_threshold: Minimum score threshold as a fraction of optimal score (0.0-1.0)
        pre_computed_matrix: Optional pre-computed score matrix to skip computation

    Returns:
        List of tuples (aligned_seq1, aligned_seq2, score) sorted by score
    """
    # Convert to uppercase for consistency
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Initialize the scoring matrix and traceback matrix
    m, n = len(seq1), len(seq2)
    
    if pre_computed_matrix is not None:
        # Use the pre-computed matrix
        score_matrix = pre_computed_matrix
        optimal_score = score_matrix[m][n]
    else:
        # Compute the score matrix
        score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
        
        # Initialize consecutive G tracking matrix
        g_count_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

        # Initialize first row and column with gap penalties
        for i in range(1, m + 1):
            score_matrix[i][0] = score_matrix[i - 1][0] - 1

        for j in range(1, n + 1):
            score_matrix[0][j] = score_matrix[0][j - 1] - 1

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
                score_matrix[i][j] = max(diagonal_score, up_score, left_score)
                
        optimal_score = score_matrix[m][n]
    
    # Store all possible moves at each cell
    traceback = [[[] for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Compute the traceback matrix based on the score matrix
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 and j > 0:
                traceback[i][j].append(("left", score_matrix[i][j]))
            elif j == 0 and i > 0:
                traceback[i][j].append(("up", score_matrix[i][j]))
            elif i > 0 and j > 0:
                # Calculate scores for each possible move
                diagonal_score = score_matrix[i-1][j-1]
                up_score = score_matrix[i-1][j]
                left_score = score_matrix[i][j-1]
                
                # Find the best score at this position
                best_score = score_matrix[i][j]
                min_acceptable = best_score * score_threshold
                
                # Add moves that are within threshold
                if i > 0 and j > 0 and diagonal_score + calculate_score(seq1, seq2, i-1, j-1, 0) >= min_acceptable:
                    traceback[i][j].append(("diagonal", diagonal_score))
                if i > 0 and up_score - 1 >= min_acceptable:
                    traceback[i][j].append(("up", up_score))
                if j > 0 and left_score - 1 >= min_acceptable:
                    traceback[i][j].append(("left", left_score))
                
                # Sort moves by score (highest first)
                traceback[i][j].sort(key=lambda x: x[1], reverse=True)

    # Generate multiple alignments using backtracking
    alignments = []
    min_score_threshold = optimal_score * score_threshold

    # Use depth-first search to find multiple paths
    def backtrack(i, j, aligned1, aligned2, path_score):
        if i == 0 and j == 0:
            # We've reached the beginning of both sequences
            alignments.append(
                (
                    "".join(reversed(aligned1)),
                    "".join(reversed(aligned2)),
                    path_score,
                )
            )
            return

        # Try all possible moves from this cell
        for move, move_score in traceback[i][j]:
            if move == "diagonal":
                # Calculate the score for this specific match/mismatch
                if i > 0 and j > 0:
                    # Check for consecutive Gs
                    consecutive_g_bonus = 0
                    if (
                        seq1[i - 1] == "G"
                        and seq2[j - 1] == "G"
                        and i > 1
                        and j > 1
                        and seq1[i - 2] == "G"
                        and seq2[j - 2] == "G"
                    ):
                        # Find how many consecutive Gs we have
                        k = 2
                        while (
                            i - k >= 0
                            and j - k >= 0
                            and seq1[i - k] == "G"
                            and seq2[j - k] == "G"
                        ):
                            k += 1
                        consecutive_g_bonus = k - 1

                    # Calculate position-specific score
                    pos_score = 2  # Basic match score
                    if seq1[i - 1] == seq2[j - 1]:  # Match
                        if seq1[i - 1] == "G":  # G match
                            pos_score += 1 + consecutive_g_bonus
                    else:  # Mismatch
                        pos_score = -1

                    backtrack(
                        i - 1,
                        j - 1,
                        aligned1 + [seq1[i - 1]],
                        aligned2 + [seq2[j - 1]],
                        path_score + pos_score,
                    )
                else:
                    backtrack(
                        i - 1,
                        j - 1,
                        aligned1 + [seq1[i - 1]],
                        aligned2 + [seq2[j - 1]],
                        path_score,
                    )
            elif move == "up":
                backtrack(
                    i - 1,
                    j,
                    aligned1 + [seq1[i - 1]],
                    aligned2 + ["-"],
                    path_score - 1,  # Gap penalty
                )
            elif move == "left":
                backtrack(
                    i,
                    j - 1,
                    aligned1 + ["-"],
                    aligned2 + [seq2[j - 1]],
                    path_score - 1,  # Gap penalty
                )

            # Stop if we have enough alignments
            if (
                len(alignments) >= num_alignments * 3
            ):  # Get more than needed, then filter
                return

    # Start backtracking from the bottom-right cell
    backtrack(m, n, [], [], 0)  # Start with score 0 and calculate during backtracking

    # Calculate actual scores for each alignment
    scored_alignments = []
    for aligned_seq1, aligned_seq2, path_score in alignments:
        # Recalculate score to ensure accuracy
        final_score = calculate_alignment_score(aligned_seq1, aligned_seq2)
        scored_alignments.append((aligned_seq1, aligned_seq2, final_score))

    # Sort alignments by score (highest first) and remove duplicates
    unique_alignments = []
    seen = set()

    # Find the optimal score after recalculation
    if scored_alignments:
        max_score = max(score for _, _, score in scored_alignments)
        min_acceptable_score = max_score * score_threshold
    else:
        min_acceptable_score = 0

    for aligned_seq1, aligned_seq2, score in sorted(
        scored_alignments, key=lambda x: x[2], reverse=True
    ):
        alignment_key = (aligned_seq1, aligned_seq2)
        if alignment_key not in seen and score >= min_acceptable_score:
            seen.add(alignment_key)
            unique_alignments.append((aligned_seq1, aligned_seq2, score))

            # Stop if we have enough unique alignments
            if len(unique_alignments) >= num_alignments:
                break

    return unique_alignments


def calculate_alignment_score(aligned_seq1, aligned_seq2):
    """
    Calculate the score for an alignment based on matches, mismatches, gaps,
    and consecutive G bonuses.

    Args:
        aligned_seq1, aligned_seq2: The aligned sequences

    Returns:
        The total alignment score
    """
    score = 0
    consecutive_g_count = 0

    for i in range(min(len(aligned_seq1), len(aligned_seq2))):
        if aligned_seq1[i] == aligned_seq2[i]:
            # Match
            base_score = 2
            if aligned_seq1[i] == "G":
                # G match bonus
                base_score += 1
                # Consecutive G bonus
                consecutive_g_count += 1
                base_score += (
                    consecutive_g_count - 1
                )  # Additional bonus for consecutive Gs
            else:
                consecutive_g_count = 0
            score += base_score
        else:
            # Mismatch
            score -= 1
            consecutive_g_count = 0

    return score


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
        if i < len(aligned_seq2) and aligned_seq1[i] == aligned_seq2[i]:
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
        print(
            f"Longest consecutive G run: {max(consecutive_g_runs) if consecutive_g_runs else 0}"
        )


def read_quadruplexes_from_directory(directory_path):
    """
    Read all JSON files in a directory and extract quadruplex objects.

    Args:
        directory_path: Path to directory containing JSON files

    Returns:
        List of tuples (quadruplex, source_file)
    """
    import os

    all_quadruplexes = []

    try:
        # Check if directory exists
        if not os.path.isdir(directory_path):
            print(f"Error: Directory '{directory_path}' not found.")
            return []

        # Get all JSON files in the directory
        json_files = [
            f for f in os.listdir(directory_path) if f.lower().endswith(".json")
        ]

        if not json_files:
            print(f"Warning: No JSON files found in directory '{directory_path}'.")
            return []

        print(f"Found {len(json_files)} JSON files in '{directory_path}'.")

        # Process each JSON file
        for json_file in json_files:
            file_path = os.path.join(directory_path, json_file)
            quadruplexes = read_quadruplex_json(file_path)

            # Add source file information to each quadruplex
            for quad in quadruplexes:
                all_quadruplexes.append((quad, json_file))

        print(f"Loaded {len(all_quadruplexes)} quadruplex structures in total.")
        return all_quadruplexes

    except Exception as e:
        print(f"Error reading directory: {str(e)}")
        return []


def display_quadruplex_details(quadruplex):
    """
    Display detailed information about a quadruplex, including segment analysis.

    Args:
        quadruplex: QuadruplexDotBracket object
    """
    print(quadruplex)

    # Display segment information
    segments = quadruplex.get_segments()
    if len(segments) > 1:
        print("\nSegment Analysis:")
        for i, (seq, struct, chi, loop) in enumerate(segments, 1):
            print(f"  Segment #{i}:")
            print(f"    Sequence: {seq} (length: {len(seq)})")
            print(f"    Structure: {struct}")
            print(f"    Chi: {chi}")
            print(f"    Loop: {loop}")


def compute_alignment_score_matrix(seq1, seq2):
    """
    Compute the alignment score matrix for two sequences.
    
    Args:
        seq1, seq2: The sequences to align
        
    Returns:
        Tuple of (score_matrix, optimal_score)
    """
    # Convert to uppercase for consistency
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Initialize the scoring matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Initialize consecutive G tracking matrix
    g_count_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # Initialize first row and column with gap penalties
    for i in range(1, m + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] - 1

    for j in range(1, n + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] - 1

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
            score_matrix[i][j] = max(diagonal_score, up_score, left_score)
    
    # Return the score matrix and the optimal score
    return score_matrix, score_matrix[m][n]


def align_against_quadruplexes(
    sequence, quadruplexes, num_alignments=5, score_threshold=0.8
):
    """
    Align a sequence against all quadruplex sequences and rank the results.

    Args:
        sequence: The sequence to align
        quadruplexes: List of tuples (quadruplex, source_file)
        num_alignments: Number of alignments to generate per quadruplex
        score_threshold: Score threshold for alignments

    Returns:
        List of tuples (quadruplex, source_file, aligned_seq1, aligned_seq2, score)
        sorted by score (highest first)
    """
    # First, compute score matrices for all quadruplexes
    print("Computing alignment scores for all quadruplexes...")
    
    quad_scores = []
    for quad, source_file in quadruplexes:
        # Compute score matrix and optimal score
        score_matrix, optimal_score = compute_alignment_score_matrix(sequence, quad.sequence)
        quad_scores.append((quad, source_file, score_matrix, optimal_score))
    
    # Find the best overall score
    if not quad_scores:
        return []
        
    best_overall_score = max(score for _, _, _, score in quad_scores)
    min_acceptable_score = best_overall_score * score_threshold
    
    print(f"Best overall alignment score: {best_overall_score}")
    print(f"Minimum acceptable score: {min_acceptable_score}")
    
    # Filter quadruplexes by score threshold
    filtered_quads = [
        (quad, source_file, score_matrix) 
        for quad, source_file, score_matrix, score in quad_scores
        if score >= min_acceptable_score
    ]
    
    print(f"Found {len(filtered_quads)} quadruplexes with scores above threshold.")
    
    # Generate alignments only for filtered quadruplexes
    all_alignments = []
    
    for quad, source_file, _ in filtered_quads:
        print(f"Generating alignments for sequence from {source_file}...")
        
        # Align the sequence against the quadruplex sequence using pre-computed matrix
        alignments = align_sequences(
            sequence, quad.sequence, num_alignments, score_threshold, _
        )
        
        # Add quadruplex and source information to each alignment
        for aligned_seq1, aligned_seq2, score in alignments:
            all_alignments.append(
                (
                    quad,
                    source_file,
                    0,  # No segment number
                    aligned_seq1,
                    aligned_seq2,
                    score,
                )
            )
    
    # Sort all alignments by score (highest first)
    all_alignments.sort(key=lambda x: x[5], reverse=True)
    
    return all_alignments


def display_all_alignments(alignments):
    """
    Display multiple alignments with their scores.

    Args:
        alignments: List of tuples (aligned_seq1, aligned_seq2, score)
    """
    print(f"\nFound {len(alignments)} alignment(s):")

    for i, (aligned_seq1, aligned_seq2, score) in enumerate(alignments, 1):
        display_alignment(aligned_seq1, aligned_seq2, score, i)


def display_ranked_alignments(ranked_alignments, top_n=10):
    """
    Display ranked alignments against quadruplexes.

    Args:
        ranked_alignments: List of tuples (quadruplex, source_file, segment_num, aligned_seq1, aligned_seq2, score)
        top_n: Number of top results to display
    """
    if not ranked_alignments:
        print("\nNo alignments found.")
        return

    # Limit to top N results
    results_to_display = ranked_alignments[:top_n]

    print(
        f"\nTop {len(results_to_display)} alignments (out of {len(ranked_alignments)} total):"
    )

    for i, (
        quad,
        source_file,
        _,  # Unused segment_num
        aligned_seq1,
        aligned_seq2,
        score,
    ) in enumerate(results_to_display, 1):
        print(f"\nRank #{i} (Score: {score}):")
        print(f"Source: {source_file}")
        print(f"Quadruplex: {quad.sequence}")

        # Display the alignment
        print(f"Sequence:   {aligned_seq1}")

        # Create a match line
        match_line = ""
        for j in range(len(aligned_seq1)):
            if j < len(aligned_seq2) and aligned_seq1[j] == aligned_seq2[j]:
                if aligned_seq1[j] == "G":
                    match_line += "*"  # Special indicator for G matches
                else:
                    match_line += "|"  # Regular match
            else:
                match_line += " "  # Mismatch or gap

        print(f"            {match_line}")
        print(f"Quadruplex: {aligned_seq2}")

        # Count G matches
        g_matches = sum(
            1
            for j in range(min(len(aligned_seq1), len(aligned_seq2)))
            if aligned_seq1[j] == "G" and aligned_seq2[j] == "G"
        )
        print(f"G matches: {g_matches}")


def main():
    """Main function to run the alignment tool."""
    args = parse_arguments()

    # Get the sequence to align
    sequence = args.sequence

    # Validate sequence
    if not validate_sequence(sequence):
        print(
            "Error: Invalid sequence. Please use only A, T, G, C, U, or N characters."
        )
        sys.exit(1)

    # Print input sequence
    print(f"Input Sequence: {sequence}")

    # Read quadruplexes from directory
    quadruplexes = read_quadruplexes_from_directory(args.directory)

    if not quadruplexes:
        print("No valid quadruplex structures found. Exiting.")
        sys.exit(1)

    # Align sequence against all quadruplexes
    print(f"Aligning sequence against {len(quadruplexes)} quadruplex structures...")
    print(
        f"Using parameters: num_alignments={args.num_alignments}, score_threshold={args.score_threshold}"
    )

    ranked_alignments = align_against_quadruplexes(
        sequence, quadruplexes, args.num_alignments, args.score_threshold
    )

    # Display top ranked alignments
    display_ranked_alignments(ranked_alignments, args.top_results)


if __name__ == "__main__":
    main()
