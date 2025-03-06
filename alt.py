#!/usr/bin/env python3
"""
Alternative DNA/RNA Quadruplex Analysis Tool

This program provides an alternative implementation for working with quadruplex structures.
"""

import argparse
import json
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass


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
        Validate that all fields have the same length and that the structure
        doesn't contain consecutive identical letters.

        Returns:
            bool: True if valid, False otherwise
        """
        # Check if all fields have the same length
        fields = [self.sequence, self.structure, self.chi, self.loop]
        if len(set(len(field) for field in fields)) != 1:
            return False

        # Check for consecutive identical letters in structure
        prev_char = None
        for char in self.structure:
            if char.isalpha() and char.lower() == prev_char:
                # Found consecutive identical letters (e.g., "qq")
                return False
            prev_char = char.lower() if char.isalpha() else None

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
                if quad:
                    if quad.validate():
                        quadruplexes.append(quad)
                    else:
                        # Check specifically for consecutive letters
                        has_consecutive = False
                        prev_char = None
                        for char in quad.structure:
                            if char.isalpha() and char.lower() == prev_char:
                                has_consecutive = True
                                break
                            prev_char = char.lower() if char.isalpha() else None

                        if has_consecutive:
                            print_warning(
                                f"Skipping quadruplex in {file_path} - structure contains consecutive identical letters"
                            )
                        else:
                            print_warning(
                                f"Invalid quadruplexDotBracket object found in {file_path} (fields have different lengths)"
                            )
                else:
                    print_warning(
                        f"Failed to parse quadruplexDotBracket object in {file_path}"
                    )
            else:
                # Try to parse as a direct quadruplex object
                quad = parse_quadruplex_object(data)
                if quad:
                    if quad.validate():
                        quadruplexes.append(quad)
                    else:
                        # Check specifically for consecutive letters
                        has_consecutive = False
                        prev_char = None
                        for char in quad.structure:
                            if char.isalpha() and char.lower() == prev_char:
                                has_consecutive = True
                                break
                            prev_char = char.lower() if char.isalpha() else None

                        if has_consecutive:
                            print_warning(
                                f"Skipping quadruplex in {file_path} - structure contains consecutive identical letters"
                            )
                        else:
                            print_warning(
                                f"Invalid quadruplex object found in {file_path} (fields have different lengths)"
                            )
                else:
                    print_warning(f"Failed to parse quadruplex object in {file_path}")
        elif isinstance(data, list):
            # Array of objects - could be array of quadruplexDotBracket containers or direct objects
            for item in data:
                if isinstance(item, dict) and "quadruplexDotBracket" in item:
                    quad_data = item["quadruplexDotBracket"]
                    quad = parse_quadruplex_object(quad_data)
                else:
                    quad = parse_quadruplex_object(item)

                if quad:
                    if quad.validate():
                        quadruplexes.append(quad)
                    else:
                        # Check specifically for consecutive letters
                        has_consecutive = False
                        prev_char = None
                        for char in quad.structure:
                            if char.isalpha() and char.lower() == prev_char:
                                has_consecutive = True
                                break
                            prev_char = char.lower() if char.isalpha() else None

                        if has_consecutive:
                            print_warning(
                                f"Skipping quadruplex in {file_path} - structure contains consecutive identical letters"
                            )
                        else:
                            print_warning(
                                f"Invalid quadruplex object found in {file_path} (fields have different lengths)"
                            )
                else:
                    print_warning(f"Failed to parse quadruplex object in {file_path}")

        return quadruplexes

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.", file=sys.stderr)
        return []
    except json.JSONDecodeError:
        print(f"Error: File '{file_path}' contains invalid JSON.", file=sys.stderr)
        return []
    except Exception as e:
        print(f"Error reading quadruplex data: {str(e)}", file=sys.stderr)
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
            print(f"Error: Directory '{directory_path}' not found.", file=sys.stderr)
            return []

        # Get all JSON files in the directory
        json_files = [
            f for f in os.listdir(directory_path) if f.lower().endswith(".json")
        ]

        if not json_files:
            print_warning(f"No JSON files found in directory '{directory_path}'.")
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
        print(f"Error reading directory: {str(e)}", file=sys.stderr)
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

                # Create compressed representation
                compressed_repr = compress_structure(str_representation)

                # Add position list, string representation, and compressed representation to the result
                all_combinations.append(
                    (display_combo, str_representation, compressed_repr)
                )

    return all_combinations


def compress_structure(structure):
    """
    Compress a structure string into a list of letter groups.
    Ignores dots and only keeps consecutive letters of the same type.

    Args:
        structure: The structure string to compress

    Returns:
        list: List of letter groups (e.g., 'qr', 'rq', etc.)
    """
    if not structure:
        return []

    # Normalize the structure first (make letters lowercase, remove linkers)
    normalized = []
    for char in structure:
        if char.isalpha():
            normalized.append(char.lower())
        elif char not in "-&":  # Skip linkers, keep other characters as dots
            normalized.append(".")

    normalized = "".join(normalized)

    # Extract letter groups
    result = []
    current_group = ""

    for char in normalized:
        if char.isalpha():
            current_group += char
        elif current_group:
            # We hit a non-letter after a letter group
            result.append(current_group)
            current_group = ""

    # Add the last group if there is one
    if current_group:
        result.append(current_group)

    return result


def concatenate_groups(groups):
    """
    Concatenate a list of letter groups into a single string.

    Args:
        groups: List of letter groups

    Returns:
        str: Concatenated string of all groups
    """
    return "".join(groups)


def get_unique_letters(groups):
    """
    Get the set of unique letters used in the groups.

    Args:
        groups: List of letter groups or a concatenated string

    Returns:
        set: Set of unique letters
    """
    if isinstance(groups, list):
        # Flatten the list of groups into a single string
        all_letters = "".join(groups)
    else:
        all_letters = groups

    return set(all_letters)


def compare_concatenated_structures(concat1, concat2):
    """
    Compare two concatenated structure strings directly.

    Args:
        concat1, concat2: Two concatenated structure strings

    Returns:
        float: Similarity score between 0.0 and 1.0
    """
    if not concat1 or not concat2:
        return 0.0

    # Calculate similarity based on longest common subsequence
    lcs_length = longest_common_subsequence(concat1, concat2)

    # Normalize by the length of the longer string
    max_length = max(len(concat1), len(concat2))
    if max_length == 0:
        return 0.0

    similarity = lcs_length / max_length
    return similarity


def try_letter_mapping(concat1, concat2, mapping):
    """
    Apply a letter mapping to concat2 and compare with concat1.

    Args:
        concat1: First concatenated structure string
        concat2: Second concatenated structure string
        mapping: Dictionary mapping letters in concat2 to letters in concat1

    Returns:
        float: Similarity score after mapping
    """
    # Apply the mapping to concat2
    mapped_concat2 = "".join(mapping.get(c, c) for c in concat2)

    # Compare the mapped string with concat1
    return compare_concatenated_structures(concat1, mapped_concat2)


def generate_letter_mappings(letters1, letters2):
    """
    Generate all possible mappings from letters2 to letters1.

    Args:
        letters1: Set of letters in the first structure
        letters2: Set of letters in the second structure

    Returns:
        list: List of mapping dictionaries
    """
    import itertools

    # If letters2 has fewer unique letters than letters1, we can't map
    if len(letters2) < len(letters1):
        return []

    # If they have the same number of letters, there's only one mapping to try
    if len(letters2) == len(letters1):
        # Try all permutations of letters1 to match with letters2
        mappings = []
        for perm in itertools.permutations(letters1):
            mapping = {l2: l1 for l2, l1 in zip(sorted(letters2), perm)}
            mappings.append(mapping)
        return mappings

    # If letters2 has more letters, we need to try different subsets
    mappings = []

    # For each possible subset of letters2 of size len(letters1)
    for subset in itertools.combinations(letters2, len(letters1)):
        # For each permutation of letters1
        for perm in itertools.permutations(letters1):
            # Create a mapping from the subset to letters1
            mapping = {l2: l1 for l2, l1 in zip(subset, perm)}

            # For remaining letters in letters2, map to the closest letter in letters1
            remaining = letters2 - set(subset)
            for l2 in remaining:
                # Map to the first letter in letters1 (arbitrary choice)
                mapping[l2] = next(iter(letters1))

            mappings.append(mapping)

    return mappings


def compare_compressed_structures(groups1, groups2):
    """
    Compare two compressed structure representations using concatenation
    and letter mapping approach.

    Args:
        groups1, groups2: Lists of letter groups

    Returns:
        float: Similarity score between 0.0 and 1.0, or -1.0 if incompatible
    """
    if not groups1 or not groups2:
        return 0.0

    # Concatenate the groups
    concat1 = concatenate_groups(groups1)
    concat2 = concatenate_groups(groups2)

    # Get unique letters in each structure
    letters1 = get_unique_letters(concat1)
    letters2 = get_unique_letters(concat2)

    # If the second structure has fewer unique letters, it can't match
    if len(letters2) < len(letters1):
        return -1.0  # Incompatible

    # Try all possible letter mappings
    mappings = generate_letter_mappings(letters1, letters2)

    # Find the best mapping
    best_score = 0.0
    for mapping in mappings:
        score = try_letter_mapping(concat1, concat2, mapping)
        best_score = max(best_score, score)

    return best_score


def longest_common_subsequence_groups(groups1, groups2):
    """
    Find the length of the longest common subsequence between two lists of groups.

    Args:
        groups1, groups2: Two lists of letter groups

    Returns:
        int: Length of the longest common subsequence
    """
    m, n = len(groups1), len(groups2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Groups are considered matching if they share at least one letter
            if has_common_letters(groups1[i - 1], groups2[j - 1]):
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    return dp[m][n]


def has_common_letters(group1, group2):
    """
    Check if two letter groups have at least one letter in common.

    Args:
        group1, group2: Two letter groups

    Returns:
        bool: True if they share at least one letter, False otherwise
    """
    return bool(set(group1) & set(group2))


def find_matching_groups(groups1, groups2):
    """
    Find pairs of matching groups between two lists.

    Args:
        groups1, groups2: Two lists of letter groups

    Returns:
        list: List of tuples (group1, group2) of matching groups
    """
    matches = []

    # Use dynamic programming to find the matching groups
    m, n = len(groups1), len(groups2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if has_common_letters(groups1[i - 1], groups2[j - 1]):
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    # Backtrack to find the matching pairs
    i, j = m, n
    while i > 0 and j > 0:
        if has_common_letters(groups1[i - 1], groups2[j - 1]):
            matches.append((groups1[i - 1], groups2[j - 1]))
            i -= 1
            j -= 1
        elif dp[i - 1][j] > dp[i][j - 1]:
            i -= 1
        else:
            j -= 1

    # Reverse to get the matches in the correct order
    matches.reverse()

    return matches


def calculate_group_similarity(group1, group2):
    """
    Calculate similarity between two letter groups.

    Args:
        group1, group2: Two letter groups

    Returns:
        float: Similarity score between 0.0 and 1.0
    """
    # Count common letters
    common_letters = set(group1) & set(group2)

    # Calculate Jaccard similarity
    union_size = len(set(group1) | set(group2))
    if union_size == 0:
        return 0.0

    jaccard = len(common_letters) / union_size

    # Consider order similarity
    order_similarity = longest_common_subsequence(group1, group2) / max(
        len(group1), len(group2)
    )

    # Combine the two measures
    return (jaccard * 0.5) + (order_similarity * 0.5)


def generate_alignment(combination_repr, quadruplex, mapping=None):
    """
    Generate an alignment-like string between the combination and quadruplex.

    Args:
        combination_repr: String representation of a tetrad combination
        quadruplex: QuadruplexDotBracket object
        mapping: Optional dictionary mapping letters in quadruplex to letters in combination

    Returns:
        tuple: (aligned_combo, aligned_quad, match_line) where:
            - aligned_combo is the combination with gaps inserted
            - aligned_quad is the quadruplex structure with gaps inserted
            - match_line shows matches between the two
    """
    # Compress both structures into letter groups
    combo_groups = compress_structure(combination_repr)
    quad_groups = compress_structure(quadruplex.structure)

    # Get the concatenated representations
    combo_concat = concatenate_groups(combo_groups)
    quad_concat = concatenate_groups(quad_groups)

    # If a mapping is provided, apply it to the quadruplex
    if mapping:
        mapped_quad_concat = "".join(mapping.get(c, c) for c in quad_concat)
    else:
        mapped_quad_concat = quad_concat

    # Create alignment strings
    aligned_combo = []
    aligned_quad = []
    match_line = []

    # Track positions in the original strings
    i, j = 0, 0

    # Find the longest common subsequence path
    lcs_path = get_lcs_path(combo_concat, mapped_quad_concat)

    for move in lcs_path:
        if move == "match":
            # Both strings advance with a match
            aligned_combo.append(combo_concat[i])
            aligned_quad.append(quad_concat[j])
            match_line.append("|")
            i += 1
            j += 1
        elif move == "combo_gap":
            # Insert a gap in the combination
            aligned_combo.append("-")
            aligned_quad.append(quad_concat[j])
            match_line.append(" ")
            j += 1
        elif move == "quad_gap":
            # Insert a gap in the quadruplex
            aligned_combo.append(combo_concat[i])
            aligned_quad.append("-")
            match_line.append(" ")
            i += 1

    # Convert to strings
    aligned_combo = "".join(aligned_combo)
    aligned_quad = "".join(aligned_quad)
    match_line = "".join(match_line)

    return aligned_combo, aligned_quad, match_line


def get_lcs_path(str1, str2):
    """
    Get the path of operations to transform str1 into str2 using LCS.

    Args:
        str1, str2: Two strings

    Returns:
        list: List of operations ('match', 'combo_gap', or 'quad_gap')
    """
    m, n = len(str1), len(str2)

    # Build the LCS matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    # Backtrack to find the path
    path = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and str1[i - 1] == str2[j - 1]:
            path.append("match")
            i -= 1
            j -= 1
        elif j > 0 and (i == 0 or dp[i][j - 1] >= dp[i - 1][j]):
            path.append("combo_gap")
            j -= 1
        elif i > 0 and (j == 0 or dp[i][j - 1] < dp[i - 1][j]):
            path.append("quad_gap")
            i -= 1

    # Reverse the path to get the correct order
    path.reverse()
    return path


def compare_combination_to_quadruplex(combination_repr, quadruplex):
    """
    Compare a tetrad combination representation against a QuadruplexDotBracket object.
    Uses concatenated compressed structure representation with letter mapping.

    Args:
        combination_repr: String representation of a tetrad combination
        quadruplex: QuadruplexDotBracket object

    Returns:
        float: Similarity score between 0.0 and 1.0, or -1.0 if incompatible
    """
    # Compress both structures into letter groups
    combo_groups = compress_structure(combination_repr)
    quad_groups = compress_structure(quadruplex.structure)

    # Get unique letters in each structure
    combo_letters = get_unique_letters(combo_groups)
    quad_letters = get_unique_letters(quad_groups)

    # If the quadruplex has fewer unique letters, it can't match
    if len(quad_letters) < len(combo_letters):
        return -1.0  # Incompatible

    # Compare the compressed structures
    return compare_compressed_structures(combo_groups, quad_groups)


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


def compare_compressed_lists(list1, list2):
    """
    Compare two compressed structure lists directly without concatenation.

    Args:
        list1, list2: Two lists of letter groups

    Returns:
        float: Similarity score between 0.0 and 1.0
    """
    if not list1 or not list2:
        return 0.0

    # Calculate similarity based on the longest common subsequence of groups
    lcs_length = longest_common_subsequence_groups(list1, list2)

    # Normalize by the length of the longer list
    max_length = max(len(list1), len(list2))
    if max_length == 0:
        return 0.0

    # Base similarity on LCS length
    base_similarity = lcs_length / max_length

    # Add bonus for matching group content
    content_similarity = 0.0
    if lcs_length > 0:
        # Find matching groups and calculate their content similarity
        matches = find_matching_groups(list1, list2)
        if matches:
            content_scores = []
            for g1, g2 in matches:
                # Score based on letter-by-letter similarity
                group_sim = calculate_group_similarity(g1, g2)
                content_scores.append(group_sim)

            # Average the content similarity scores
            content_similarity = sum(content_scores) / len(content_scores)

    # Combine base similarity with content similarity
    final_score = (base_similarity * 0.6) + (content_similarity * 0.4)

    return min(1.0, final_score)  # Cap at 1.0


def process_combination_comparison(args):
    """
    Process a single combination comparison against multiple quadruplexes.
    This function is designed to be used with multiprocessing.

    Args:
        args: Tuple containing (combination_index, str_repr, compressed_repr, combo, quadruplexes)

    Returns:
        Tuple of (combination_index, combo, str_repr, compressed_repr, list of (similarity, compressed_similarity, quad_index, quad, sources, alignment))
    """
    combination_index, str_repr, compressed_repr, combo, quadruplexes = args

    # Compare this combination to each quadruplex
    quad_scores = []
    for quad_index, (quad, sources) in enumerate(quadruplexes):
        # First filter by concatenated similarity
        concat_similarity = compare_combination_to_quadruplex(str_repr, quad)

        # Only include perfect matches (score of 1.0) for concatenated comparison
        if concat_similarity == 1.0:
            # Calculate compressed list similarity for ranking
            quad_compressed = compress_structure(quad.structure)
            compressed_similarity = compare_compressed_lists(
                compressed_repr, quad_compressed
            )

            # Get unique letters in each structure
            combo_concat = concatenate_groups(compressed_repr)
            quad_concat = concatenate_groups(quad_compressed)
            combo_letters = get_unique_letters(combo_concat)
            quad_letters = get_unique_letters(quad_concat)

            # Find the best mapping for alignment
            best_mapping = None
            best_score = -1

            if len(quad_letters) >= len(combo_letters):
                mappings = generate_letter_mappings(combo_letters, quad_letters)
                for mapping in mappings:
                    score = try_letter_mapping(combo_concat, quad_concat, mapping)
                    if score > best_score:
                        best_score = score
                        best_mapping = mapping

            # Generate alignment
            alignment = generate_alignment(str_repr, quad, best_mapping)

            quad_scores.append(
                (
                    concat_similarity,
                    compressed_similarity,
                    quad_index,
                    quad,
                    sources,
                    alignment,
                )
            )

    # Sort by compressed similarity score (highest first)
    quad_scores.sort(key=lambda x: x[1], reverse=True)

    return (combination_index, combo, str_repr, compressed_repr, quad_scores)


def generate_sequence_alignment(
    sequence, quadruplex, combination_repr, structure_alignment
):
    """
    Generate a detailed alignment between the original sequence and quadruplex sequence.

    Args:
        sequence: Original DNA/RNA sequence
        quadruplex: QuadruplexDotBracket object
        combination_repr: String representation of the tetrad combination
        structure_alignment: Tuple of (aligned_combo, aligned_quad, match_line)

    Returns:
        tuple: (aligned_seq, aligned_quad_seq, match_line) showing the sequence alignment
    """
    aligned_combo, aligned_quad, match_line = structure_alignment

    # Create a mapping from positions in combination_repr to positions in sequence
    combo_to_seq = {}
    seq_pos = 0
    for i, char in enumerate(combination_repr):
        if char.isalpha():  # This is a G position
            combo_to_seq[i] = seq_pos
            seq_pos += 1
        elif char == ".":  # This is a non-G position
            combo_to_seq[i] = seq_pos
            seq_pos += 1

    # Create a mapping from positions in quad.structure to positions in quad.sequence
    quad_to_seq = {}
    seq_pos = 0
    for i, char in enumerate(quadruplex.structure):
        if char != "-" and char != "&":  # Skip linkers
            quad_to_seq[i] = seq_pos
            seq_pos += 1

    # Generate the sequence alignment
    aligned_seq = []
    aligned_quad_seq = []
    seq_match_line = []

    combo_pos = 0
    quad_pos = 0

    for i in range(len(aligned_combo)):
        if aligned_combo[i] == "-":
            # Gap in combination
            aligned_seq.append("-")
            if aligned_quad[i] == "-":
                # Gap in both (shouldn't happen)
                aligned_quad_seq.append("-")
                seq_match_line.append(" ")
            else:
                # Get the corresponding position in quadruplex sequence
                struct_pos = quad_pos
                if struct_pos < len(quadruplex.structure):
                    seq_pos = quad_to_seq.get(struct_pos, None)
                    if seq_pos is not None and seq_pos < len(quadruplex.sequence):
                        aligned_quad_seq.append(quadruplex.sequence[seq_pos])
                    else:
                        aligned_quad_seq.append("?")
                else:
                    aligned_quad_seq.append("?")
                seq_match_line.append(" ")
                quad_pos += 1
        elif aligned_quad[i] == "-":
            # Gap in quadruplex
            if combo_pos < len(combination_repr):
                seq_pos = combo_to_seq.get(combo_pos, None)
                if seq_pos is not None and seq_pos < len(sequence):
                    aligned_seq.append(sequence[seq_pos])
                else:
                    aligned_seq.append("?")
            else:
                aligned_seq.append("?")
            aligned_quad_seq.append("-")
            seq_match_line.append(" ")
            combo_pos += 1
        else:
            # Match or mismatch
            if combo_pos < len(combination_repr):
                seq_pos = combo_to_seq.get(combo_pos, None)
                if seq_pos is not None and seq_pos < len(sequence):
                    aligned_seq.append(sequence[seq_pos])
                else:
                    aligned_seq.append("?")
            else:
                aligned_seq.append("?")

            if quad_pos < len(quadruplex.structure):
                seq_pos = quad_to_seq.get(quad_pos, None)
                if seq_pos is not None and seq_pos < len(quadruplex.sequence):
                    aligned_quad_seq.append(quadruplex.sequence[seq_pos])
                else:
                    aligned_quad_seq.append("?")
            else:
                aligned_quad_seq.append("?")

            # Check if the nucleotides match
            if aligned_seq[-1] == aligned_quad_seq[-1]:
                seq_match_line.append("|")
            elif aligned_seq[-1] in "TU" and aligned_quad_seq[-1] in "TU":
                seq_match_line.append("|")  # T and U are considered matches
            else:
                seq_match_line.append(" ")

            combo_pos += 1
            quad_pos += 1

    return "".join(aligned_seq), "".join(aligned_quad_seq), "".join(seq_match_line)


def find_best_matches_parallel(
    tetrad_combinations,
    quadruplexes,
):
    """
    Find the best matches between tetrad combinations and quadruplexes in parallel.
    Processes all combinations against all quadruplexes.

    Args:
        tetrad_combinations: List of tetrad combinations
        quadruplexes: List of quadruplexes
        num_results: Number of top results to return

    Returns:
        List of (combination_index, combo, str_repr, list of top matches)
    """
    if not tetrad_combinations or not quadruplexes:
        return []

    # Process all combinations against all quadruplexes
    # Prepare arguments for parallel processing
    process_args = [
        (i, str_repr, compressed_repr, combo, quadruplexes)
        for i, (combo, str_repr, compressed_repr) in enumerate(tetrad_combinations)
    ]

    # Determine the number of processes to use
    num_processes = min(multiprocessing.cpu_count(), len(process_args))
    if num_processes < 1:
        num_processes = 1

    print(
        f"Processing {len(tetrad_combinations)} combinations against {len(quadruplexes)} quadruplexes using {num_processes} processes..."
    )

    # Process comparisons in parallel
    results = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        for result in executor.map(process_combination_comparison, process_args):
            combination_index, combo, str_repr, compressed_repr, quad_scores = result
            # Keep all matches
            if quad_scores:  # Only add combinations that have matches
                results.append(
                    (combination_index, combo, str_repr, compressed_repr, quad_scores)
                )
            # Show progress
            print(
                f"  Processed combination {combination_index + 1}/{len(tetrad_combinations)}",
                end="\r",
            )

    print()  # New line after progress

    # Sort results by the best compressed similarity score (highest first)
    results.sort(key=lambda x: x[4][0][1] if x[4] and x[4][0] else 0, reverse=True)

    return results


def print_warning(message):
    """Print a warning message to stderr."""
    print(f"Warning: {message}", file=sys.stderr)


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
        "-r",
        "--results",
        type=int,
        default=10,
        help="Number of top results to display in final ranking (default: 10)",
    )
    parser.add_argument(
        "--show-combinations",
        action="store_true",
        help="Display all generated tetrad combinations",
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
            "Error: Invalid sequence. Please use only A, T, G, C, U, or N characters.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Print input sequence
    seq_type = "RNA" if is_rna else "DNA"
    print(f"Input Sequence: {sequence} ({seq_type})")

    # Calculate and display maximum possible tetrads
    max_tetrads = calculate_max_tetrads(sequence)
    print(f"Maximum possible tetrads (G-quartets): {max_tetrads}")

    # Generate tetrad combinations
    all_tetrad_combinations = generate_tetrad_combinations(sequence)
    print(f"\nFound {len(all_tetrad_combinations)} possible tetrad combinations")

    # Filter to unique compressed representations
    unique_compressed = {}
    for combo, str_repr, compressed_repr in all_tetrad_combinations:
        # Convert list to tuple for hashability
        compressed_key = tuple(compressed_repr)
        if compressed_key not in unique_compressed:
            unique_compressed[compressed_key] = (combo, str_repr, compressed_repr)

    # Convert back to list
    tetrad_combinations = list(unique_compressed.values())
    print(
        f"After filtering duplicates: {len(tetrad_combinations)} unique tetrad patterns"
    )

    # Display the combinations if requested
    if args.show_combinations:
        display_limit = min(10, len(tetrad_combinations))  # Default to showing 10
        for i, (combo, str_repr, compressed_repr) in enumerate(
            tetrad_combinations[:display_limit], 1
        ):
            print(f"\nCombination {i}:")
            print(f"  Representation: {str_repr}")
            print(f"  Compressed: {compressed_repr}")
            print(f"  Concatenated: {concatenate_groups(compressed_repr)}")
            print(f"  Unique letters: {sorted(get_unique_letters(compressed_repr))}")
            for j, tetrad in enumerate(combo, 1):
                print(f"  Tetrad {j}: G positions {', '.join(map(str, tetrad))}")

        if len(tetrad_combinations) > display_limit:
            print(
                f"\n... and {len(tetrad_combinations) - display_limit} more combinations"
            )

    # Read quadruplexes from directory
    quadruplexes = read_quadruplexes_from_directory(args.directory)

    if not quadruplexes:
        print("No valid quadruplex structures found. Exiting.", file=sys.stderr)
        sys.exit(1)

    # For now, just print summary of loaded quadruplexes
    dna_count = sum(1 for quad, _ in quadruplexes if not quad.is_rna)
    rna_count = sum(1 for quad, _ in quadruplexes if quad.is_rna)

    print("\nSummary of loaded quadruplexes:")
    print(f"  Total unique structures: {len(quadruplexes)}")
    print(f"  DNA structures: {dna_count}")
    print(f"  RNA structures: {rna_count}")

    # If we have both tetrad combinations and quadruplexes, compare them
    if tetrad_combinations and quadruplexes:
        print("\nComparing tetrad combinations to quadruplex structures...")

        # Find best matches in parallel
        best_matches = find_best_matches_parallel(
            tetrad_combinations,
            quadruplexes,
        )

        # Sort results by the best compressed similarity score (highest first)
        best_matches.sort(
            key=lambda x: x[4][0][1] if x[4] and x[4][0] else 0, reverse=True
        )

        # Limit to requested number of results
        display_results = best_matches[: min(args.results, len(best_matches))]

        # Display results
        print(
            f"\nTop {len(display_results)} combinations with perfect matches (out of {len(best_matches)} total):"
        )

        if not display_results:
            print("  No perfect matches found.")

        for rank, (
            combination_index,
            combo,
            str_repr,
            compressed_repr,
            quad_scores,
        ) in enumerate(display_results, 1):
            print(f"\nRank #{rank} - Combination {combination_index + 1}: {str_repr}")
            print(f"  Compressed: {compressed_repr}")

            # Display tetrad positions
            for j, tetrad in enumerate(combo, 1):
                print(f"  Tetrad {j}: G positions {', '.join(map(str, tetrad))}")

            # Display top matches
            print("  Top matches:")
            combo_concat = concatenate_groups(compressed_repr)
            combo_letters = get_unique_letters(combo_concat)

            for (
                concat_similarity,
                compressed_similarity,
                quad_index,
                quad,
                sources,
                alignment,
            ) in quad_scores:
                quad_compressed = compress_structure(quad.structure)
                quad_concat = concatenate_groups(quad_compressed)
                quad_letters = get_unique_letters(quad_concat)

                print(f"    Concat match score: {concat_similarity:.2f}")
                print(f"    Compressed match score: {compressed_similarity:.2f}")
                print(f"    Structure: {quad.structure}")
                print(f"    Compressed: {quad_compressed}")
                print(f"    Concatenated: {quad_concat}")
                print(f"    Unique letters: {sorted(quad_letters)}")

                # If the quadruplex has more letters, show the best mapping
                if len(quad_letters) > len(combo_letters):
                    mappings = generate_letter_mappings(combo_letters, quad_letters)
                    best_mapping = None
                    best_score = -1

                    for mapping in mappings:
                        score = try_letter_mapping(combo_concat, quad_concat, mapping)
                        if score > best_score:
                            best_score = score
                            best_mapping = mapping

                    if best_mapping:
                        print(f"    Best mapping: {best_mapping}")
                        mapped_concat = "".join(
                            best_mapping.get(c, c) for c in quad_concat
                        )
                        print(f"    Mapped: {mapped_concat}")

                # Display the alignment
                aligned_combo, aligned_quad, match_line = alignment
                print(f"    Structure Alignment:")
                print(f"      Combination: {aligned_combo}")
                print(f"      Match:       {match_line}")
                print(f"      Quadruplex:  {aligned_quad}")

                # Generate and display sequence alignment
                seq_alignment = generate_sequence_alignment(
                    sequence, quad, str_repr, alignment
                )
                aligned_seq, aligned_quad_seq, seq_match_line = seq_alignment

                print(f"    Sequence Alignment:")
                print(f"      Input:      {aligned_seq}")
                print(f"      Match:      {seq_match_line}")
                print(f"      Quadruplex: {aligned_quad_seq}")

                print(
                    f"    Source: {sources[0]}"
                    + (f" (and {len(sources) - 1} more)" if len(sources) > 1 else "")
                )


if __name__ == "__main__":
    main()
