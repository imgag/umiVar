# coding=utf-8

# Import Libraries
import argparse
import sys
from collections import Counter


class Homopolymer:

    # Constructor
    def __init__(self, reference_fasta_file, min_homopolymer_length, max_outliers):
        self.reference_fasta_file = reference_fasta_file
        self.min_homopolymer_length = min_homopolymer_length
        self.max_outliers = max_outliers
        self.chromosome = ''
        self.chr_sequence = ''
        self.homopolymer_positions = dict()

    # Parse the fasta file
    def parse_fasta(self):

        with open(self.reference_fasta_file) as f1:
            for line in f1:
                line = line.rstrip('\n')

                # Read Fasta header line
                if line.startswith('>'):

                    # Print homopolymers for previous chromosome
                    for key, value in sorted(self.homopolymer_positions.items(), key=lambda x: x[0]):
                        if value == 1:
                            print("{}\t{}".format(self.chromosome, key))

                    # Reset containers
                    self.chr_sequence = ''
                    self.homopolymer_positions = dict()

                    # Read new chromosome name
                    split_header_line = line.split(" ")
                    self.chromosome = split_header_line[0][1:]
                    self.chr_sequence = ''

                else:
                    self.chr_sequence += line

        self.find_homopolymers()

        # Print homopolymers of the last chromosome
        for key, value in sorted(self.homopolymer_positions.items(), key=lambda x: x[0]):
            if value == 1:
                print("{}\t{}".format(self.chromosome, key))

    # Find homopolymer sequences in a chromosome
    def find_homopolymers(self):
        position = 1
        current_seq = ''

        # Parse sequence of chromosome
        for i in self.chr_sequence:

            # Store last n nucleotides of reference sequence
            if len(current_seq) < self.min_homopolymer_length:
                current_seq += i
            else:
                # Check if last n bases are a homopolymer
                homopolymer = self.count_max_nucleotide(current_seq)

                # If homopolymer, label all positions overlapping homopolymer in dictionary
                if homopolymer == 1:
                    for pos in range(position - self.min_homopolymer_length, position):
                        self.homopolymer_positions[pos] = homopolymer

                # Shift sequence by one nucleotide
                current_seq = current_seq[1:]
                current_seq += i

            position += 1

    # Identify the most common nucleotide in sequence and label homopolymers
    def count_max_nucleotide(self, sequence):
        wc = Counter(sequence)
        most_common_base = wc.most_common(1)[0][1]

        # Check if most common nucleotide forms a homopolymer
        if most_common_base >= (self.min_homopolymer_length - self.max_outliers):
            return 1
        else:
            return 0


# Run homopolymer finder
def main():

    # Read parameters
    parser = argparse.ArgumentParser(description='Find all homopolymers of length n in reference sequence (fasta).')
    parser.add_argument('-r', '--ref', type=str, help='Reference fasta file.', required=True)
    parser.add_argument('-l', '--length', type=int, required=False, default=6,
                        help='Minimum length of homopolymers to report. Default = 6')
    parser.add_argument('-o', '--outlier', type=int, required=False, default=0,
                        help='Maximum number of outliers (other nucleotides) in homopolymers. Default = 0')

    try:
        args = parser.parse_args()
    except IOError as io:
        print(io)
        sys.exit('Error reading parameters.')

    reference_fasta_file = args.ref
    min_homopolymer_length = args.length
    max_outliers = args.outlier

    # Create Homopolymer object and run parse_fasta method
    parser = Homopolymer(reference_fasta_file, min_homopolymer_length, max_outliers)
    parser.parse_fasta()


# Run homopolymer finder
main()
