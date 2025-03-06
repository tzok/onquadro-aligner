#!/usr/bin/env python3
"""
Alternative DNA/RNA Quadruplex Analysis Tool

This program provides an alternative implementation for working with quadruplex structures.
"""

import sys
import argparse
import json
import os
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


if __name__ == "__main__":
    main()
