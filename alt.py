#!/usr/bin/env python3
"""
Alternative DNA/RNA Quadruplex Analysis Tool

This program provides an alternative implementation for working with quadruplex structures.
"""

import sys
import argparse
import json
import os
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
    is_rna: bool = False

    def __str__(self):
        """String representation of the quadruplex."""
        seq_type = "RNA" if self.is_rna else "DNA"
        return (
            f"Sequence: {self.sequence} ({seq_type})\n"
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

        # Determine if this is RNA or DNA
        is_rna = "U" in sequence and "T" not in sequence

        # Create and return the object
        return QuadruplexDotBracket(sequence, structure, chi, loop, is_rna)

    except Exception as e:
        print(f"Error parsing quadruplex object: {str(e)}")
        return None


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


def read_quadruplexes_from_directory(directory_path):
    """
    Read all JSON files in a directory and extract quadruplex objects.
    Groups duplicate quadruplexes (same sequence and structure) together.

    Args:
        directory_path: Path to directory containing JSON files

    Returns:
        List of tuples (quadruplex, source_files) where source_files is a list of file names
    """
    # Dictionary to store unique quadruplexes by their sequence and structure
    unique_quadruplexes = {}  # Key: (sequence, structure), Value: (quad, [source_files])

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
        total_quads_loaded = 0
        for json_file in json_files:
            file_path = os.path.join(directory_path, json_file)
            quadruplexes = read_quadruplex_json(file_path)
            total_quads_loaded += len(quadruplexes)

            # Group quadruplexes by sequence and structure
            for quad in quadruplexes:
                key = (quad.sequence, quad.structure)
                if key in unique_quadruplexes:
                    # Add this file to the existing quadruplex's sources
                    unique_quadruplexes[key][1].append(json_file)
                else:
                    # Create a new entry
                    unique_quadruplexes[key] = (quad, [json_file])

        # Convert dictionary to list of tuples
        all_quadruplexes = [
            (quad, source_files)
            for (quad, source_files) in unique_quadruplexes.values()
        ]

        print(f"Loaded {total_quads_loaded} quadruplex structures in total.")
        print(
            f"After grouping duplicates: {len(all_quadruplexes)} unique quadruplex structures."
        )
        return all_quadruplexes

    except Exception as e:
        print(f"Error reading directory: {str(e)}")
        return []


def calculate_max_tetrads(sequence):
    """
    Calculate the theoretical maximum number of tetrads possible in a sequence.

    Args:
        sequence: The DNA/RNA sequence

    Returns:
        int: Maximum number of tetrads (G-quartets) possible
    """
    # Count the number of G nucleotides in the sequence
    g_count = sequence.upper().count("G")

    # Each tetrad requires 4 G nucleotides
    max_tetrads = g_count // 4

    return max_tetrads


def find_g_positions(sequence):
    """
    Find the positions of all G nucleotides in a sequence.

    Args:
        sequence: The DNA/RNA sequence

    Returns:
        list: List of positions (0-based) where G nucleotides are found
    """
    return [i for i, nucleotide in enumerate(sequence.upper()) if nucleotide == "G"]


def generate_tetrad_combinations(sequence):
    """
    Generate all possible combinations of tetrad assignments for a sequence.

    Rules:
    1. Two consecutive Gs cannot be part of the same tetrad
    2. Each tetrad contains exactly four Gs

    Args:
        sequence: The DNA/RNA sequence

    Returns:
        list: List of valid tetrad combinations, where each combination is a tuple containing:
              - a list of tetrads (each tetrad is a tuple of four G positions, 1-based)
              - a string representation with letters at tetrad positions and dots elsewhere
    """
    # Find positions of all G nucleotides
    g_positions = find_g_positions(sequence)

    # If we have fewer than 4 Gs, no tetrads are possible
    if len(g_positions) < 4:
        return []

    # Generate all possible combinations of 4 G positions
    from itertools import combinations

    # Store valid tetrads
    valid_tetrads = []

    # Check each possible combination of 4 G positions
    for combo in combinations(g_positions, 4):
        # Check if any consecutive Gs are in this combination
        is_valid = True
        for i in range(len(g_positions) - 1):
            if (
                g_positions[i] in combo
                and g_positions[i + 1] in combo
                and g_positions[i + 1] == g_positions[i] + 1
            ):
                is_valid = False
                break

        if is_valid:
            valid_tetrads.append(combo)

    # If we have fewer than 2 valid tetrads, no complete combinations are possible
    if len(valid_tetrads) < 1:
        return []

    # Generate all possible combinations of tetrads
    from itertools import combinations as itercombo

    # Maximum number of tetrads possible
    max_tetrads = len(g_positions) // 4

    all_combinations = []

    # Try different numbers of tetrads, from max down to 1
    for num_tetrads in range(max_tetrads, 0, -1):
        for tetrad_combo in itercombo(valid_tetrads, num_tetrads):
            # Check if each G is used at most once
            used_positions = set()
            is_valid = True

            for tetrad in tetrad_combo:
                for pos in tetrad:
                    if pos in used_positions:
                        is_valid = False
                        break
                    used_positions.add(pos)

                if not is_valid:
                    break

            if is_valid:
                # Convert to 1-based positions for display
                display_combo = [
                    tuple(pos + 1 for pos in tetrad) for tetrad in tetrad_combo
                ]

                # Create string representation
                # Start with all dots
                str_representation = ["." for _ in range(len(sequence))]

                # Fill in tetrad positions with appropriate letters
                for i, tetrad in enumerate(tetrad_combo):
                    # Get the letter for this tetrad (q, r, s, t, ...)
                    tetrad_letter = chr(ord("q") + i)

                    # Place the letter at each G position in this tetrad
                    for pos in tetrad:
                        # Convert from 0-based to string index
                        str_representation[pos] = tetrad_letter

                # Join the list into a string
                str_representation = "".join(str_representation)

                # Add both the position list and string representation to the result
                all_combinations.append((display_combo, str_representation))

    return all_combinations


def normalize_structure(structure):
    """
    Normalize a structure string by making all letters lowercase,
    removing linkers, and replacing non-letters with dots.

    Args:
        structure: The structure string to normalize

    Returns:
        str: Normalized structure string
    """
    result = []
    for char in structure:
        if char.isalpha():
            result.append(char.lower())
        elif char in "-&":
            # Skip linkers
            continue
        else:
            result.append(".")

    return "".join(result)


def compare_combination_to_quadruplex(combination_repr, quadruplex):
    """
    Compare a tetrad combination representation against a QuadruplexDotBracket object.

    Args:
        combination_repr: String representation of a tetrad combination
        quadruplex: QuadruplexDotBracket object

    Returns:
        float: Similarity score between 0.0 and 1.0
    """
    # Normalize the quadruplex structure
    norm_structure = normalize_structure(quadruplex.structure)

    # If lengths don't match, we need to align them
    if len(combination_repr) != len(norm_structure):
        return align_and_score_structures(combination_repr, norm_structure)

    # If lengths match, we can do a direct comparison
    return score_matching_structures(combination_repr, norm_structure)


def score_matching_structures(str1, str2):
    """
    Score two structures of the same length.

    Args:
        str1, str2: Two structure strings of the same length

    Returns:
        float: Similarity score between 0.0 and 1.0
    """
    if len(str1) != len(str2):
        raise ValueError("Structures must be of the same length")

    total_positions = len(str1)
    if total_positions == 0:
        return 0.0

    matches = 0
    letter_positions = []

    # First pass: identify letter positions and exact matches
    for i in range(total_positions):
        # Both are letters - check if they match
        if str1[i].isalpha() and str2[i].isalpha():
            if str1[i] == str2[i]:
                matches += 1.0  # Full match
            else:
                matches += 0.2  # Different letters, but both are letters
            letter_positions.append(i)
        # One is a letter, one is a dot
        elif str1[i].isalpha() or str2[i].isalpha():
            matches += 0.1  # Small partial match for letter vs dot
        # Both are dots
        else:
            matches += 0.5  # Dots match dots with medium weight

    # Weight letter positions more heavily
    letter_weight = 2.0 if letter_positions else 1.0

    # Calculate final score
    score = matches / (total_positions * letter_weight)
    return min(1.0, score)  # Cap at 1.0


def align_and_score_structures(str1, str2):
    """
    Align two structures of different lengths and score their similarity.
    Uses a simplified alignment approach that prioritizes matching letter positions.

    Args:
        str1, str2: Two structure strings that may have different lengths

    Returns:
        float: Similarity score between 0.0 and 1.0
    """
    # Extract letter patterns (ignore dots)
    pattern1 = "".join([c for c in str1 if c.isalpha()])
    pattern2 = "".join([c for c in str2 if c.isalpha()])

    if not pattern1 or not pattern2:
        # If either has no letters, score based on proportion of dots
        dot_ratio = min(str1.count(".") / len(str1), str2.count(".") / len(str2))
        return dot_ratio * 0.5  # Lower score for dot-only matches

    # Find longest common subsequence of letters
    lcs_length = longest_common_subsequence(pattern1, pattern2)

    # Calculate letter pattern similarity
    pattern_similarity = lcs_length / max(len(pattern1), len(pattern2))

    # Calculate overall structure similarity
    # Consider both letter pattern and overall length
    length_ratio = min(len(str1), len(str2)) / max(len(str1), len(str2))

    # Weight pattern similarity more heavily
    final_score = (pattern_similarity * 0.8) + (length_ratio * 0.2)

    return final_score


def longest_common_subsequence(str1, str2):
    """
    Find the length of the longest common subsequence between two strings.

    Args:
        str1, str2: Two strings

    Returns:
        int: Length of the longest common subsequence
    """
    m, n = len(str1), len(str2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    return dp[m][n]


def process_combination_comparison(args):
    """
    Process a single combination comparison against multiple quadruplexes.
    This function is designed to be used with multiprocessing.

    Args:
        args: Tuple containing (combination_index, str_repr, combo, quadruplexes)

    Returns:
        Tuple of (combination_index, combo, str_repr, list of (similarity, quad_index, quad, sources))
    """
    combination_index, str_repr, combo, quadruplexes = args
    
    # Compare this combination to each quadruplex
    quad_scores = []
    for quad_index, (quad, sources) in enumerate(quadruplexes):
        similarity = compare_combination_to_quadruplex(str_repr, quad)
        quad_scores.append((similarity, quad_index, quad, sources))
    
    # Sort by similarity score (highest first)
    quad_scores.sort(reverse=True)
    
    return (combination_index, combo, str_repr, quad_scores)


def find_best_matches_parallel(tetrad_combinations, quadruplexes, num_combinations=5, num_quadruplexes=None, top_matches=3, num_processes=0):
    """
    Find the best matches between tetrad combinations and quadruplexes in parallel.
    
    Args:
        tetrad_combinations: List of tetrad combinations
        quadruplexes: List of quadruplexes
        num_combinations: Number of combinations to process
        num_quadruplexes: Number of quadruplexes to compare against (None for all)
        top_matches: Number of top matches to return for each combination
        num_processes: Number of processes to use (0 for auto-detect)
        
    Returns:
        List of (combination_index, combo, str_repr, list of top matches)
    """
    if not tetrad_combinations or not quadruplexes:
        return []
    
    # Limit the number of combinations to process
    combinations_to_process = min(num_combinations, len(tetrad_combinations))
    
    # Limit the number of quadruplexes to compare against
    if num_quadruplexes is not None:
        quadruplexes_to_compare = quadruplexes[:min(num_quadruplexes, len(quadruplexes))]
    else:
        quadruplexes_to_compare = quadruplexes
    
    # Prepare arguments for parallel processing
    process_args = [
        (i, str_repr, combo, quadruplexes_to_compare)
        for i, (combo, str_repr) in enumerate(tetrad_combinations[:combinations_to_process])
    ]
    
    # Determine the number of processes to use
    if num_processes <= 0:
        num_processes = min(multiprocessing.cpu_count(), len(process_args))
    if num_processes < 1:
        num_processes = 1
    
    print(f"Processing {combinations_to_process} combinations against {len(quadruplexes_to_compare)} quadruplexes using {num_processes} processes...")
    
    # Process comparisons in parallel
    results = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        for result in executor.map(process_combination_comparison, process_args):
            combination_index, combo, str_repr, quad_scores = result
            # Keep only the top matches
            top_quad_scores = quad_scores[:top_matches]
            results.append((combination_index, combo, str_repr, top_quad_scores))
            # Show progress
            print(f"  Processed combination {combination_index + 1}/{combinations_to_process}", end="\r")
    
    print()  # New line after progress
    
    # Sort results by the best match score (highest first)
    results.sort(key=lambda x: x[3][0][0] if x[3] else 0, reverse=True)
    
    return results


def validate_sequence(sequence):
    """
    Validate if the input is a valid DNA/RNA sequence.

    Args:
        sequence: The sequence to validate

    Returns:
        Tuple of (is_valid, is_rna) where is_valid is a boolean indicating if the sequence is valid,
        and is_rna is a boolean indicating if the sequence is RNA (True) or DNA (False)
    """
    valid_chars = set("ATGCUN-")  # Include '-' as a valid character
    if not all(char.upper() in valid_chars for char in sequence):
        return False, False

    # Determine if this is RNA or DNA
    is_rna = "U" in sequence and "T" not in sequence

    return True, is_rna


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Alternative tool for working with quadruplex structures."
    )
    parser.add_argument("sequence", help="DNA/RNA sequence to analyze")
    parser.add_argument(
        "-d",
        "--directory",
        help="Directory containing JSON files with quadruplex data",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--combinations",
        type=int,
        default=10,
        help="Number of tetrad combinations to process (default: 10)",
    )
    parser.add_argument(
        "-q",
        "--quadruplexes",
        type=int,
        default=5,
        help="Number of quadruplexes to compare against (default: 5)",
    )
    parser.add_argument(
        "--all-quadruplexes",
        action="store_true",
        help="Compare against all quadruplexes instead of just the first N",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=0,
        help="Number of processes to use (default: auto-detect)",
    )
    return parser.parse_args()


def main():
    """Main function to run the tool."""
    args = parse_arguments()

    # Get the sequence to analyze
    sequence = args.sequence

    # Validate sequence
    is_valid, is_rna = validate_sequence(sequence)
    if not is_valid:
        print(
            "Error: Invalid sequence. Please use only A, T, G, C, U, or N characters."
        )
        sys.exit(1)

    # Print input sequence
    seq_type = "RNA" if is_rna else "DNA"
    print(f"Input Sequence: {sequence} ({seq_type})")

    # Calculate and display maximum possible tetrads
    max_tetrads = calculate_max_tetrads(sequence)
    print(f"Maximum possible tetrads (G-quartets): {max_tetrads}")

    # Generate and display tetrad combinations
    tetrad_combinations = generate_tetrad_combinations(sequence)
    print(f"\nFound {len(tetrad_combinations)} possible tetrad combinations:")

    # Display the combinations (limit based on user preference)
    display_limit = min(args.combinations, len(tetrad_combinations))
    for i, (combo, str_repr) in enumerate(tetrad_combinations[:display_limit], 1):
        print(f"\nCombination {i}:")
        print(f"  Representation: {str_repr}")
        for j, tetrad in enumerate(combo, 1):
            print(f"  Tetrad {j}: G positions {', '.join(map(str, tetrad))}")

    if len(tetrad_combinations) > display_limit:
        print(f"\n... and {len(tetrad_combinations) - display_limit} more combinations")

    # Read quadruplexes from directory
    quadruplexes = read_quadruplexes_from_directory(args.directory)

    if not quadruplexes:
        print("No valid quadruplex structures found. Exiting.")
        sys.exit(1)

    # For now, just print summary of loaded quadruplexes
    dna_count = sum(1 for quad, _ in quadruplexes if not quad.is_rna)
    rna_count = sum(1 for quad, _ in quadruplexes if quad.is_rna)

    print(f"\nSummary of loaded quadruplexes:")
    print(f"  Total unique structures: {len(quadruplexes)}")
    print(f"  DNA structures: {dna_count}")
    print(f"  RNA structures: {rna_count}")

    # If we have both tetrad combinations and quadruplexes, compare them
    if tetrad_combinations and quadruplexes:
        print("\nComparing tetrad combinations to quadruplex structures...")
        
        # Set multiprocessing start method if on macOS
        if sys.platform == "darwin":  # macOS
            multiprocessing.set_start_method("spawn", force=True)
        
        # Determine which quadruplexes to use
        quad_limit = None if args.all_quadruplexes else args.quadruplexes
        
        # Find best matches in parallel
        best_matches = find_best_matches_parallel(
            tetrad_combinations, 
            quadruplexes,
            num_combinations=args.combinations,
            num_quadruplexes=quad_limit,
            top_matches=3,
            num_processes=args.processes
        )
        
        # Display results
        print(f"\nTop {len(best_matches)} combinations with best matches:")
        
        for rank, (combination_index, combo, str_repr, quad_scores) in enumerate(best_matches, 1):
            print(f"\nRank #{rank} - Combination {combination_index + 1}: {str_repr}")
            
            # Display tetrad positions
            for j, tetrad in enumerate(combo, 1):
                print(f"  Tetrad {j}: G positions {', '.join(map(str, tetrad))}")
            
            # Display top matches
            print("  Top matches:")
            for similarity, quad_index, quad, sources in quad_scores:
                print(f"    Match score: {similarity:.2f}")
                print(f"    Structure: {quad.structure}")
                print(f"    Normalized: {normalize_structure(quad.structure)}")
                print(
                    f"    Source: {sources[0]}"
                    + (f" (and {len(sources) - 1} more)" if len(sources) > 1 else "")
                )


if __name__ == "__main__":
    main()
