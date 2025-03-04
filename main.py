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


def align_sequences(seq1, seq2):
    """
    Align two DNA/RNA sequences.

    This is a placeholder function that will be implemented in future steps.
    Currently, it just returns the original sequences.
    """
    # Placeholder for alignment algorithm
    # Will be implemented in future steps
    aligned_seq1 = seq1
    aligned_seq2 = seq2

    return aligned_seq1, aligned_seq2


def display_alignment(aligned_seq1, aligned_seq2):
    """Display the aligned sequences."""
    print("\nAlignment Result:")
    print(f"Sequence 1: {aligned_seq1}")
    print(f"Sequence 2: {aligned_seq2}")


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
