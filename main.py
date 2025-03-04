#!/usr/bin/env python3
"""
DNA/RNA Sequence Alignment Tool

This program aligns two DNA or RNA sequences provided as command line arguments.
"""

import sys
import argparse
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import List, Optional, Tuple


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
        Get the segments of the quadruplex (parts separated by '&').

        Returns:
            List of tuples (sequence_segment, structure_segment, chi_segment, loop_segment)
        """
        seq_segments = self.sequence.split("&")
        struct_segments = self.structure.split("&")
        chi_segments = self.chi.split("&")
        loop_segments = self.loop.split("&")

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
        sequence = str(data.get("sequence", "")).replace("-", "&")
        structure = str(data.get("structure", "")).replace("-", "&")

        # Handle chi field - ensure it's a string
        chi_data = data.get("chi", "")
        if not isinstance(chi_data, str):
            chi = str(chi_data).replace("-", "&")
        else:
            chi = chi_data.replace("-", "&")

        # Handle loop field - ensure it's a string
        loop_data = data.get("loop", "")
        if not isinstance(loop_data, str):
            if isinstance(loop_data, list):
                loop = "&".join(str(x) for x in loop_data)
            else:
                loop = str(loop_data).replace("-", "&")
        else:
            loop = loop_data.replace("-", "&")

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


def align_sequences(seq1, seq2, score_threshold=0.8, pre_computed_matrix=None):
    """
    Align two DNA/RNA sequences with emphasis on consecutive G matches.
    Treats T and U as matches.

    Uses dynamic programming with beam search to find multiple optimal and suboptimal
    alignments that maximize the score, with special consideration for G nucleotides.

    Args:
        seq1, seq2: The sequences to align
        score_threshold: Minimum score threshold as a fraction of optimal score (0.0-1.0)
        pre_computed_matrix: Optional pre-computed score matrix to skip computation

    Returns:
        List of tuples (aligned_seq1, aligned_seq2, score) sorted by score
    """
    # Convert to uppercase for consistency
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Initialize the scoring matrix
    m, n = len(seq1), len(seq2)

    if pre_computed_matrix is not None:
        # Use the pre-computed matrix
        score_matrix = pre_computed_matrix
        optimal_score = score_matrix[m][n]
    else:
        # Compute the score matrix
        score_matrix, optimal_score = compute_alignment_score_matrix(seq1, seq2)

    # Calculate the minimum acceptable score
    min_score_threshold = optimal_score * score_threshold

    # Initialize a list to store alignments
    alignments = []

    # Use beam search to find multiple alignments
    # Each state is (i, j, aligned_seq1, aligned_seq2, current_score, g_count)
    # where i, j are positions in the matrix, aligned_seq1/2 are the alignments so far,
    # current_score is the score so far, and g_count is the number of consecutive G matches

    # Start with the bottom-right cell
    beam = [(m, n, [], [], 0, 0)]  # (i, j, aligned1, aligned2, score, g_count)
    seen_states = set()  # To avoid revisiting states

    # Keep track of complete alignments
    complete_alignments = []

    # Beam search parameters
    beam_width = 50  # Keep this many paths at each step
    max_iterations = m * n * 2  # Limit total iterations to avoid excessive computation

    iteration = 0
    while beam and iteration < max_iterations:
        iteration += 1

        # Get the current state
        i, j, aligned1, aligned2, current_score, g_count = beam.pop(0)

        # If we've reached the beginning of both sequences, we have a complete alignment
        if i == 0 and j == 0:
            # Convert lists to strings
            aligned_seq1 = "".join(aligned1)
            aligned_seq2 = "".join(aligned2)

            # Add to complete alignments
            complete_alignments.append((aligned_seq1, aligned_seq2, current_score))

            # Continue to process other paths
            continue

        # Generate possible moves from this cell
        moves = []

        # Diagonal move (match/mismatch)
        if i > 0 and j > 0:
            # Check for G match and consecutive Gs
            is_g_match = seq1[i - 1] == "G" and seq2[j - 1] == "G"
            consecutive_g_bonus = 0
            new_g_count = 0

            if is_g_match:
                new_g_count = g_count + 1
                consecutive_g_bonus = g_count  # Bonus based on previous consecutive Gs

            # Calculate score for this move
            # Check for exact match or T-U match
            is_match = seq1[i - 1] == seq2[j - 1] or (
                seq1[i - 1] in "TU" and seq2[j - 1] in "TU"
            )

            if is_match:  # Match
                move_score = 2  # Basic match score
                if is_g_match:
                    move_score += 1 + consecutive_g_bonus  # G match bonus
            else:  # Mismatch
                move_score = -1

            # Add diagonal move
            moves.append(
                (
                    i - 1,
                    j - 1,
                    [seq1[i - 1]] + aligned1,
                    [seq2[j - 1]] + aligned2,
                    current_score + move_score,
                    new_g_count,
                )
            )

        # Up move (gap in seq2)
        if i > 0:
            moves.append(
                (
                    i - 1,
                    j,
                    [seq1[i - 1]] + aligned1,
                    ["-"] + aligned2,
                    current_score - 1,  # Gap penalty
                    0,  # Reset consecutive G count
                )
            )

        # Left move (gap in seq1)
        if j > 0:
            moves.append(
                (
                    i,
                    j - 1,
                    ["-"] + aligned1,
                    [seq2[j - 1]] + aligned2,
                    current_score - 1,  # Gap penalty
                    0,  # Reset consecutive G count
                )
            )

        # Add valid moves to the beam
        for move in moves:
            # Create a state key that doesn't include the alignments (to save memory)
            state_key = (move[0], move[1])

            # Only add if we haven't seen this state before or if the score is better
            if state_key not in seen_states:
                beam.append(move)
                seen_states.add(state_key)

        # Sort the beam by score (highest first) and keep only the top paths
        beam.sort(key=lambda x: x[4], reverse=True)
        beam = beam[:beam_width]

    # Calculate actual scores for each alignment
    scored_alignments = []
    for aligned_seq1, aligned_seq2, _ in complete_alignments:
        # Calculate the final score
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

    return unique_alignments


def calculate_alignment_score(aligned_seq1, aligned_seq2):
    """
    Calculate the score for an alignment based on matches, mismatches, gaps,
    and consecutive G bonuses. Treats T and U as matches.

    Args:
        aligned_seq1, aligned_seq2: The aligned sequences

    Returns:
        The total alignment score
    """
    score = 0
    consecutive_g_count = 0

    for i in range(min(len(aligned_seq1), len(aligned_seq2))):
        # Check for exact match or T-U match
        is_match = aligned_seq1[i] == aligned_seq2[i] or (
            aligned_seq1[i] in "TU" and aligned_seq2[i] in "TU"
        )

        if is_match:
            # Match
            base_score = 2
            if aligned_seq1[i] == "G" and aligned_seq2[i] == "G":
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
    Treats T and U as matches.

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
        if i < len(aligned_seq2):
            # Check for exact match or T-U match
            is_match = aligned_seq1[i] == aligned_seq2[i] or (
                aligned_seq1[i] in "TU" and aligned_seq2[i] in "TU"
            )

            if is_match:
                if aligned_seq1[i] == "G" and aligned_seq2[i] == "G":
                    match_line += "*"  # Special indicator for G matches
                else:
                    match_line += "|"  # Regular match
            else:
                match_line += " "  # Mismatch or gap
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


def process_quadruplex(args):
    """
    Process a single quadruplex for alignment score computation.
    This function is designed to be used with multiprocessing.

    Args:
        args: Tuple containing (sequence, quad, source_file)

    Returns:
        Tuple of (quad, source_file, score_matrix, optimal_score)
    """
    sequence, quad, source_file = args
    score_matrix, optimal_score = compute_alignment_score_matrix(
        sequence, quad.sequence
    )
    return (quad, source_file, score_matrix, optimal_score)


def process_alignment(args):
    """
    Process a single quadruplex for alignment generation.
    This function is designed to be used with multiprocessing.

    Args:
        args: Tuple containing (sequence, quad, source_file, score_matrix, score_threshold)

    Returns:
        List of tuples (quad, source_file, 0, aligned_seq1, aligned_seq2, score)
    """
    sequence, quad, source_file, score_matrix, score_threshold = args

    # Align the sequence against the quadruplex sequence using pre-computed matrix
    alignments = align_sequences(sequence, quad.sequence, score_threshold, score_matrix)

    # Add quadruplex and source information to each alignment
    return [
        (quad, source_file, 0, aligned_seq1, aligned_seq2, score)
        for aligned_seq1, aligned_seq2, score in alignments
    ]


def compute_alignment_score_matrix(seq1, seq2):
    """
    Compute the alignment score matrix for two sequences.
    Treats T and U as matches.

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
                # Check for exact match or T-U match
                is_match = seq1[i - 1] == seq2[j - 1] or (
                    seq1[i - 1] in "TU" and seq2[j - 1] in "TU"
                )

                if is_match:
                    diagonal_score = score_matrix[i - 1][j - 1] + calculate_score(
                        seq1, seq2, i - 1, j - 1, consecutive_g_count
                    )
                else:
                    diagonal_score = score_matrix[i - 1][j - 1] - 1  # Mismatch penalty
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


def align_against_quadruplexes(sequence, quadruplexes, score_threshold=0.8):
    """
    Align a sequence against all quadruplex sequences and rank the results.
    Uses parallel processing to speed up computation.

    Args:
        sequence: The sequence to align
        quadruplexes: List of tuples (quadruplex, source_file)
        score_threshold: Score threshold for alignments

    Returns:
        List of tuples (quadruplex, source_file, aligned_seq1, aligned_seq2, score)
        sorted by score (highest first)
    """
    # Determine the number of processes to use
    num_processes = min(multiprocessing.cpu_count(), len(quadruplexes))
    if num_processes < 1:
        num_processes = 1

    print(
        f"Computing alignment scores for all quadruplexes using {num_processes} processes..."
    )

    # Prepare arguments for parallel processing
    process_args = [(sequence, quad, source_file) for quad, source_file in quadruplexes]

    # Compute score matrices in parallel
    quad_scores = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        for i, result in enumerate(executor.map(process_quadruplex, process_args)):
            quad_scores.append(result)
            # Show progress
            if i % 10 == 0 and i > 0:
                print(f"  Processed {i}/{len(quadruplexes)} quadruplexes...")

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

    # Generate alignments in parallel for filtered quadruplexes
    all_alignments = []
    total_quads = len(filtered_quads)

    if total_quads > 0:
        print(f"Generating alignments in parallel using {num_processes} processes...")

        # Prepare arguments for parallel processing
        alignment_args = [
            (sequence, quad, source_file, score_matrix, score_threshold)
            for quad, source_file, score_matrix in filtered_quads
        ]

        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            for i, result in enumerate(executor.map(process_alignment, alignment_args)):
                all_alignments.extend(result)
                # Show progress
                print(f"  Processed alignments for quadruplex {i + 1}/{total_quads}")

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
    Treats T and U as matches.

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
        # Count G matches
        g_matches = sum(
            1
            for j in range(min(len(aligned_seq1), len(aligned_seq2)))
            if aligned_seq1[j] == "G" and aligned_seq2[j] == "G"
        )

        print(f"\nRank #{i} (Score: {score}, G matches: {g_matches}):")
        print(f"Source: {source_file}")
        print(f"Quadruplex: {quad.sequence}")  # Keep ampersands

        # Display the alignment
        print(f"Sequence:   {aligned_seq1}")

        # Create a match line
        match_line = ""
        for j in range(len(aligned_seq1)):
            if j < len(aligned_seq2):
                # Check for exact match or T-U match
                is_match = aligned_seq1[j] == aligned_seq2[j] or (
                    aligned_seq1[j] in "TU" and aligned_seq2[j] in "TU"
                )

                if is_match:
                    if aligned_seq1[j] == "G" and aligned_seq2[j] == "G":
                        match_line += "*"  # Special indicator for G matches
                    else:
                        match_line += "|"  # Regular match
                else:
                    match_line += " "  # Mismatch or gap
            else:
                match_line += " "  # Mismatch or gap

        print(f"            {match_line}")
        print(f"Quadruplex: {aligned_seq2}")

        # Align structure, chi and loop with the quadruplex sequence by adding gaps
        aligned_structure = ""
        aligned_chi = ""
        aligned_loop = ""

        # Track position in the original quadruplex sequence
        orig_pos = 0

        # For each character in the aligned quadruplex sequence
        for char in aligned_seq2:
            if char == "-":
                # If it's a gap, add a gap to structure, chi and loop
                aligned_structure += "-"
                aligned_chi += "-"
                aligned_loop += "-"
            else:
                # If it's not a gap, add the corresponding character from structure, chi and loop
                if orig_pos < len(quad.structure):
                    aligned_structure += quad.structure[orig_pos]
                    aligned_chi += quad.chi[orig_pos]
                    aligned_loop += quad.loop[orig_pos]
                    orig_pos += 1
                else:
                    # Handle case where aligned sequence is longer than original
                    aligned_structure += "?"
                    aligned_chi += "?"
                    aligned_loop += "?"

        print(f"Structure:  {aligned_structure}")
        print(f"Chi:        {aligned_chi}")
        print(f"Loop:       {aligned_loop}")


def main():
    """Main function to run the alignment tool."""
    # Set multiprocessing start method
    if sys.platform == "darwin":  # macOS
        multiprocessing.set_start_method("spawn", force=True)

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
    print(f"Using parameters: score_threshold={args.score_threshold}")

    ranked_alignments = align_against_quadruplexes(
        sequence, quadruplexes, args.score_threshold
    )

    # Display top ranked alignments
    display_ranked_alignments(ranked_alignments, args.top_results)


if __name__ == "__main__":
    main()
