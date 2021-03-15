# coding=utf-8

# Import Libraries
import argparse
import sys


# Run pileup_parser
def main():
    # Read parameters
    parser = argparse.ArgumentParser(description='Find all homopolymers of length n in reference sequence (fasta).')
    parser.add_argument('-r', '--ref', type=str, help='Reference fasta file.', required=True)
    parser.add_argument('-l', '--length', type=int, required=False, default=6,
                        help='Minimum length of homopolymers to report. Default = 6')

    try:
        args = parser.parse_args()
    except IOError as io:
        print(io)
        sys.exit('Error reading parameters.')

    reference_fasta_file = args.ref
    min_homopolymer_length = args.length

    print(reference_fasta_file + "   " + str(min_homopolymer_length))


# Run homopolymer finder
main()
