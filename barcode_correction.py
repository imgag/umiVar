# coding=utf-8

# Import Libraries
import argparse
import sys
import timeit
import networkx as nx
import pysam


# Edit distance between two sequences
def similarity(a, b):
    return sum(x != y for x, y in zip(a, b))


# Return key of dictionary entry with maximum value
def key_with_max_value(d):
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


# Returns ?
def update_tag(tag, value):
    return [value if x[0] == value[0] else x for x in tag]


# Returns ?
def get_edge(code, in_list, errors):
    edge = [[code, in_list[x]] for x in range(0, len(in_list)) if (similarity(in_list[x], code) <= errors)]
    return edge


# Returns the barcode of a read
def extract_barcode(read_entry, barcode_type):
    qname = str(read_entry.qname)
    barcode_entry = qname.split(':')[-1]
    bc1 = barcode_entry.split(',')[0]
    bc2 = barcode_entry.split(',')[1]
    barcode_sequence = ""

    # Grouping by barcodes
    if barcode_type == "BEGINNING":
        barcode_sequence = bc1
    elif barcode_type == "END":
        barcode_sequence = bc2
    elif barcode_type == "BOTH":
        barcode_sequence = bc1 + bc2

    return barcode_sequence


# Returns group of reads with the same barcode
def extract_bc_groups(barcode_entries, bc_network):  # Input, list of reads with equal start and end
    barcode_groups = {}

    sorted_key_list = sorted(barcode_entries.keys(), key=lambda s: len(barcode_entries.get(s)), reverse=True)

    while len(barcode_entries) > 0:

        # The reads are stored in dict barcode_entries. Their keys are the barcodes.
        # Get the most frequent key (Barcode)
        most_frequent_barcode = sorted_key_list[0]

        # Create a new key of the most common barcode to which reads with the same barcode are added (barcode read group).
        barcode_groups[most_frequent_barcode] = list()

        sim = list(bc_network.adj[most_frequent_barcode])

        for i in sim:
            barcode_groups[most_frequent_barcode].extend(
                barcode_entries[i])  # Grouping reads based on similarity of the barcodes
            del barcode_entries[i]  # Removing barcodes already considered
            bc_network.remove_node(i)
            sorted_key_list.remove(i)

    return barcode_groups  # Dictionary with Barcode as a key and reads as values


# Return reverse complement of sequence
def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])


# Reduce mapping quality to < 20
def reduce_mapq(read_original):
    if read_original.mapq >= 20:
        read_original.mapq = 19
        read_new = read_original
    else:
        read_new = read_original

    return read_new


# Return integer value of ascii quality character
def ascii2int(ascii_values):
    quality = [ord(i) - 33 for i in ascii_values]
    return quality


# Return most frequent base at a locus
def most_common_base(nucleotides, qualities, min_bq):
    hq = [x for x in range(0, len(qualities)) if qualities[x] >= min_bq]
    result_list = [nucleotides[i] for i in hq]

    max_base = max(sorted(set(result_list)), key=result_list.count)
    num_diff_bases = len(set(result_list))

    diff_count = len([x for x in range(0, len(result_list)) if result_list[x] != max_base])
    base_count = len([x for x in range(0, len(result_list)) if result_list[x] == max_base])

    return max_base, num_diff_bases, base_count, diff_count


# Returns most frequent base of a locus disregarding quality values
def most_common_base_low_qual(nucleotides):
    max_base = max(sorted(set(nucleotides)), key=nucleotides.count)
    num_diff_bases = len(set(nucleotides))

    return max_base, num_diff_bases


# Returns qualities for a given base as integer
def get_qualities(base_of_interest, nucleotides, qualities):
    quality = [(ord(qualities[x]) - 33) for x in range(0, len(nucleotides)) if nucleotides[x] == base_of_interest]
    return quality


# Change all base qualities of the read to 0 ("!")
def error_read_qual(read_entry):
    adjusted_read = read_entry
    qualities = adjusted_read.qual
    qualities = "!" * len(qualities)
    adjusted_read.qual = qualities
    return adjusted_read


# Change discordant bases in a barcode group to Ns
def error_read_seq(read_entry):
    adjusted_read = read_entry
    adjusted_read.seq = "N" * len(adjusted_read.seq)
    return adjusted_read


# Check if there is at least one high quality base
def low_quality_base_check(qualities, min_bq):
    return max(qualities) >= min_bq


# Get consensus read qualities
def consensus_quality(qualities, min_bq, errors):

    # Get number of bases with  good quality
    copies = len([x for x in range(0, len(qualities)) if qualities[x] >= min_bq])

    # Compute base qualities of consensus base
    max_qual = max(qualities)
    if max_qual < min_bq:
        new_qual = max_qual

    else:
        # Unreliable consensus labeled with quality 0: if 3 or more discordant bases, or more than 25% discordant bases
        if (errors >= 1 and float(errors) / (copies + errors) > 0.25) or errors >= 3:
            new_qual = 0

        # Otherwise, the maximum quality value across PCR copies is used as new quality, with minimum of 30
        else:
            max_quality = max(qualities)
            new_qual = max(30, max_quality)

    return new_qual


# Compute consensus sequence of reads with the same barcode
def generate_consensus_read(reads, min_bq, set_n):
    consensus_seq = list()
    consensus_qual = list()
    consensus_read = reads[0]

    # Objects to save info in lexicographic order of the reads
    read_names_list = list()
    reads_dict = {}
    flag_dict = {}

    color = ['230,242,255', '179,215,255', '128,187,255', '77,160,255', '26,133,255']

    if len(reads) == 1:

        # Add info about the amount of duplicates per barcode family group
        count = 1
        current_color = color[0]

        # Adding barcodes to tag in bam file
        consensus_read.tags += [('DP', count)]
        consensus_read.tags += [('YC', current_color)]

        # Info about barcode groups
        log_info = (consensus_read.qname, str(consensus_read.pos), str(len(reads)))
        log_info = "\t".join(log_info) + "\n"

    else:
        # Use the characteristics of the first read's alignment to check if other reads in the barcode group
        # diverge in the number and alignment of gaps (these will be flagged as bad quality)
        first_ref_length = reads[0].reference_length
        first_read_length = reads[0].rlen
        first_cigar = reads[0].cigarstring

        # Containers for consensus read
        seq_dict = {}
        qual_dict = {}
        mapq_list = list()
        lq_read_count = 0
        last_read = reads[0]

        # Compare the sequence of reads in a barcode group
        for i in reads:

            # In case that the amount of indels differs between duplicates, take the first read as consensus, but change all base qualities to 0
            if (first_ref_length != i.reference_length) or (first_read_length != i.rlen) or (first_cigar != i.cigarstring):
                read_name = i.qname
                flag = i.flag
                reads_dict[read_name] = i
                flag_dict[read_name] = flag
                read_names_list.append(read_name)
                lq_read_count = lq_read_count + 1
                continue

            # Saving the amount of duplicates from first round of correction
            last_read = i

            # Adding reads to a dictionary to sort them in lexicographical order
            read_name = i.qname
            flag = i.flag
            reads_dict[read_name] = i
            flag_dict[read_name] = flag
            read_names_list.append(read_name)

            read_length = i.rlen
            seq = i.seq
            qual = i.qual
            mapq_list.append(i.mapq)

            soft_clip = i.pos - i.qstart

            for b in range(0, read_length):
                base = seq[b]
                base_qual = qual[b]

                real_b = soft_clip + b

                if real_b in seq_dict:
                    seq_dict[real_b].append(base)
                    qual_dict[real_b].append(base_qual)

                else:
                    seq_dict[real_b] = list()
                    qual_dict[real_b] = list()
                    seq_dict[real_b].append(base)
                    qual_dict[real_b].append(base_qual)

        for position in sorted(seq_dict):

            current_qual = ascii2int(qual_dict[position])

            if low_quality_base_check(current_qual, min_bq):
                base = most_common_base(seq_dict[position], current_qual, min_bq)

                consensus_base = base[0]
                num_diff_bases = base[1]
                consensus_base_count = base[2]
                diff_count = base[3]

                qualities = get_qualities(consensus_base, seq_dict[position], qual_dict[position])

                if num_diff_bases < 3 and consensus_base_count > diff_count:
                    consensus_quality_num = consensus_quality(qualities, min_bq, diff_count)

                    if set_n and consensus_quality_num == 0:
                        consensus_quality_num = qualities[0]
                        consensus_base = "N"

                    consensus_quality_ascii = chr(consensus_quality_num + 33)
                    consensus_qual.append(consensus_quality_ascii)
                    consensus_seq.append(consensus_base)

                elif num_diff_bases >= 3 or consensus_base_count <= diff_count:
                    consensus_quality_num = 0

                    if set_n:
                        consensus_quality_num = qualities[0]
                        consensus_base = "N"

                    consensus_quality_ascii = chr(consensus_quality_num + 33)
                    consensus_qual.append(consensus_quality_ascii)
                    consensus_seq.append(consensus_base)

                else:
                    print("Error")

            else:
                consensus_base = most_common_base_low_qual(seq_dict[position])[0]
                consensus_quality_num = 0

                if set_n:
                    qualities = get_qualities(consensus_base, seq_dict[position], qual_dict[position])
                    consensus_quality_num = qualities[0]
                    consensus_base = "N"

                consensus_seq.append(consensus_base)
                consensus_quality_ascii = chr(consensus_quality_num + 33)
                consensus_qual.append(consensus_quality_ascii)

        # Take the info from the last read in the group
        sorted_read_names = sorted(read_names_list)

        # Take as template the last HQ read, but change the read name and the flag
        consensus_read = last_read
        consensus_read.qname = sorted_read_names[0]

        # Compute average mapping quality
        if len(mapq_list) > 0:
            consensus_read.mapq = int(round(float(sum(mapq_list)) / len(mapq_list)))
        else:
            consensus_read.mapq = 0

        consensus_read.flag = flag_dict[sorted_read_names[0]]

        # Consensus seq per position
        consensus_seq = ''.join(consensus_seq)
        consensus_read.seq = consensus_seq

        # Base qualities are calculated as the mean of the base qualities of each read.
        # In case there was more than one divergent base call at the position, the consensus base quality is set to 0
        consensus_qual = ''.join(consensus_qual)
        consensus_read.qual = consensus_qual

        # Add info about the amount of duplicates per family group
        count = len(reads) - lq_read_count

        if count > 5:
            current_color = color[4]
        else:
            current_color = color[count - 1]

        # Add DP tag
        consensus_read.tags += [('DP', count)]

        # Add color
        consensus_read.tags += [('YC', current_color)]

        # Info about barcode groups
        log_info = (consensus_read.qname, str(consensus_read.pos), str(len(reads)))
        log_info = "\t".join(log_info) + "\n"

    return consensus_read, log_info


# Main method bundles argument parsing and BAM file parsing
# Example command line: python barcode_correction.py --infile PATH/TO/test.bam --outfile PATH/TO/corrected.test.bam --barcodes BOTH
def main():

    # Read script parameters
    parser = argparse.ArgumentParser(description='Correcting BAM files using barcodes info')
    parser.add_argument('--infile', required=True, dest='infile', help='Input BAM file.')
    parser.add_argument('--outfile', required=True, dest='outfile', help='Output BAM file.')
    parser.add_argument('--barcodes', required=False, dest='barcodes', choices=['START', 'END', 'BOTH'], default='BOTH',
                        help='Barcode position: START = 5\' barcode; END = 3\' barcode; BOTH = 5\' and 3\' barcodes. Default = BOTH')
    parser.add_argument('--minBQ', required=False, dest='minBQ', type=int, default=10,
                        help='Minimum base quality to be considered. Default = 30')
    parser.add_argument('--barcode_error', required=False, dest='barcode_error', type=int, default=0,
                        help='Maximum number of sequencing errors allowed in barcode sequence. Default = 0')
    parser.add_argument('--n', required=False, dest='n', action='store_true',
                        help='Use Ns instead of reducing base quality.')

    try:
        args = parser.parse_args()
    except IOError as io:
        print(io)
        sys.exit('Error reading parameters.')

    # Input BAM
    samfile = ''
    try:
        samfile = pysam.Samfile(args.infile, "rb")
    except IOError as io:
        exit("Cannot open input file. Error:\n" + io)

    # Output BAM
    outfile = ''
    try:
        outfile = pysam.Samfile(args.outfile, mode="wb", template=samfile)
    except IOError as io:
        exit("Cannot open output file. Error:\n" + io)

    logfile = open(args.outfile + ".log", 'w')

    min_bq = args.minBQ
    errors = args.barcode_error
    set_n = args.n

    pos = 0
    positions_dict = {}
    unique_barcodes = {}

    start = timeit.default_timer()

    # Parse BAM file
    for read in samfile.fetch():
        if not read.is_secondary:

            ref_start = str(read.pos)

            # Both are required. Start of next read, and tlen shows the sign of of the read (- or +), which helps to separate pair reads when they map to the same coordinates
            ref_length = str(read.next_reference_start) + ',' + str(read.tlen)

            # Getting the barcodes
            bc = extract_barcode(read, args.barcodes)  # Extract the barcode
            code = bc

            if ref_start == pos:

                # To store the codes for each ref_length
                if ref_length in positions_dict:
                    if code in positions_dict[ref_length]:
                        positions_dict[ref_length][code].append(read)
                    else:
                        positions_dict[ref_length][code] = list()
                        positions_dict[ref_length][code].append(read)
                else:
                    positions_dict[ref_length] = {}
                    positions_dict[ref_length][code] = list()
                    positions_dict[ref_length][code].append(read)

                # Allowing errors
                if errors > 0:
                    if ref_length in unique_barcodes:
                        if code in list(unique_barcodes[ref_length].nodes()):
                            unique_barcodes[ref_length].add_node(code)
                        else:
                            unique_barcodes[ref_length].add_node(code)
                            edge = get_edge(code, list(unique_barcodes[ref_length].nodes()), errors)
                            unique_barcodes[ref_length].add_edges_from(edge)
                    else:
                        unique_barcodes[ref_length] = nx.Graph()
                        unique_barcodes[ref_length].add_node(code)
                        edge = get_edge(code, list(unique_barcodes[ref_length].nodes()), errors)
                        unique_barcodes[ref_length].add_edges_from(edge)

            else:

                if len(positions_dict) > 0 and errors > 0:
                    for pos2 in positions_dict:
                        # When we allow errors in the Barcodes, we re-group them by similarity (Errors specified in parameter)
                        barcode_dict = extract_bc_groups(positions_dict[pos2], unique_barcodes[pos2])

                        for barcode in barcode_dict:
                            # Printing consensus reads to a new bam file
                            new_read, log_string = generate_consensus_read(list(barcode_dict[barcode]), min_bq, set_n)
                            logfile.write(log_string)
                            outfile.write(new_read)

                    positions_dict = {}
                    unique_barcodes = {}
                    pos = ref_start

                    if ref_length in positions_dict:
                        if code in positions_dict[ref_length]:
                            positions_dict[ref_length][code].append(read)
                        else:
                            positions_dict[ref_length][code] = list()
                            positions_dict[ref_length][code].append(read)
                    else:
                        positions_dict[ref_length] = {}
                        positions_dict[ref_length][code] = list()
                        positions_dict[ref_length][code].append(read)

                    # Allowing errors
                    if errors > 0:
                        if ref_length in unique_barcodes:
                            if code in list(unique_barcodes[ref_length].nodes()):
                                unique_barcodes[ref_length].add_node(code)
                            else:
                                unique_barcodes[ref_length].add_node(code)
                                edge = get_edge(code, list(unique_barcodes[ref_length].nodes()), errors)
                                unique_barcodes[ref_length].add_edges_from(edge)
                        else:
                            unique_barcodes[ref_length] = nx.Graph()
                            unique_barcodes[ref_length].add_node(code)
                            edge = get_edge(code, list(unique_barcodes[ref_length].nodes()), errors)
                            unique_barcodes[ref_length].add_edges_from(edge)

                elif len(positions_dict) > 0 and errors == 0:

                    barcode_dict = positions_dict

                    for pos2 in barcode_dict:
                        # printing consensus reads to a new bam file
                        for barcode in barcode_dict[pos2]:
                            new_read, log_string = generate_consensus_read(list(barcode_dict[pos2][barcode]), min_bq, set_n)

                            logfile.write(log_string)
                            outfile.write(new_read)

                    positions_dict = {}
                    pos = ref_start

                    if ref_length in positions_dict:
                        if code in positions_dict[ref_length]:
                            positions_dict[ref_length][code].append(read)
                        else:
                            positions_dict[ref_length][code] = list()
                            positions_dict[ref_length][code].append(read)
                    else:
                        positions_dict[ref_length] = {}
                        positions_dict[ref_length][code] = list()
                        positions_dict[ref_length][code].append(read)

                else:
                    positions_dict = {}
                    unique_barcodes = {}
                    pos = ref_start

                    if ref_length in positions_dict:
                        if code in positions_dict[ref_length]:
                            positions_dict[ref_length][code].append(read)
                        else:
                            positions_dict[ref_length][code] = list()
                            positions_dict[ref_length][code].append(read)
                    else:
                        positions_dict[ref_length] = {}
                        positions_dict[ref_length][code] = list()
                        positions_dict[ref_length][code].append(read)

                    # Allowing errors
                    if errors > 0:
                        if ref_length in unique_barcodes:
                            if code in list(unique_barcodes[ref_length].nodes()):
                                unique_barcodes[ref_length].add_node(code)
                            else:
                                unique_barcodes[ref_length].add_node(code)
                                edge = get_edge(code, list(unique_barcodes[ref_length].nodes()), errors)
                                unique_barcodes[ref_length].add_edges_from(edge)
                        else:
                            unique_barcodes[ref_length] = nx.Graph()
                            unique_barcodes[ref_length].add_node(code)
                            edge = get_edge(code, list(unique_barcodes[ref_length].nodes()), errors)
                            unique_barcodes[ref_length].add_edges_from(edge)

    # We need to print the last groups of reads
    if len(positions_dict) > 0 and errors > 0:
        for pos2 in positions_dict:

            # When we allow errors in the Barcodes, we re-group them by similarity (Errors specified in parameter)
            barcode_dict = extract_bc_groups(positions_dict[pos2], unique_barcodes[pos2])

            for barcode in barcode_dict:
                # Printing consensus reads to a new bam file
                new_read, log_string = generate_consensus_read(list(barcode_dict[barcode]), min_bq, set_n)

                logfile.write(log_string)
                outfile.write(new_read)

    elif len(positions_dict) > 0 and errors == 0:

        barcode_dict = positions_dict
        for pos2 in barcode_dict:

            # printing consensus reads to a new bam file
            for barcode in barcode_dict[pos2]:
                new_read, log_string = generate_consensus_read(list(barcode_dict[pos2][barcode]), min_bq, set_n)

                logfile.write(log_string)
                outfile.write(new_read)

    samfile.close()
    logfile.close()
    outfile.close()

    stop = timeit.default_timer()
    print('TIME')
    print(stop - start)


# Run program
main()
