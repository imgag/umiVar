# coding=utf-8

# Import Libraries
import argparse
import sys


# Return coverage and high quality coverage of a locus
def coverage(quality_sequence, min_bq):
    cov = len(quality_sequence)
    hq = len([x for x in range(0, len(quality_sequence)) if (ord(quality_sequence[x]) - 33) >= min_bq])
    return cov, hq


# Return strand of an observed base
def fwd_rev_strand(a):
    strand = "forward" if any(map(str.isupper, a)) else "reverse"
    return strand


# Return key of dictionary entry with maximum value
def key_with_max_value(d):
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


# Parse pileup line
def pileup_info(line, min_bq):
    line = line.rstrip('\n')

    locus, pos, ref_base, cov, bases, qualities = line.split("\t")
    pos = int(pos)
    cov = int(cov)

    # Save last base in case an indel follows (this base and the attached quality value belong to the indel)
    last_base = ''
    last_base_quality = 0

    # Reference base and N counts
    ref_count = 0
    ref_count_fwd = 0
    ref_count_rev = 0
    ref_count_lq = 0
    n_count = 0

    # Count alternative alleles
    mm_count_lq = 0
    alt_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0}
    alt_quals = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0}

    # Count indels
    indel_counts = {'del': 0, 'ins': 0}
    del_quals = {}
    ins_quals = {}

    insertion = dict()
    insertion_fwd = dict()
    insertion_rev = dict()

    deletion = dict()
    deletion_fwd = dict()
    deletion_rev = dict()

    # Loop control variables for parsing sequence and quality strings of a pileup file
    b_pos = 0
    q_pos = 0

    # Parse pileup line
    while b_pos < len(bases):
        base = bases[b_pos]

        # Match to reference
        if base in '.,':
            base_quality = ord(qualities[q_pos]) - 33
            last_base = base
            last_base_quality = base_quality

            if base_quality >= min_bq:
                ref_count += 1

                if base == ".":
                    ref_count_fwd += 1
                    alt_quals[ref_base.upper()] += base_quality

                elif base == ",":
                    ref_count_rev += 1
                    alt_quals[ref_base.lower()] += base_quality

            else:
                ref_count_lq += 1

            b_pos += 1
            q_pos += 1

        # Mismatch to reference
        elif base in 'ATCGatcg':
            base_quality = ord(qualities[q_pos]) - 33
            last_base = base
            last_base_quality = base_quality

            if base_quality >= min_bq:
                alt_counts[base] += 1
                alt_quals[base] += base_quality

            else:
                mm_count_lq += 1

            b_pos += 1
            q_pos += 1

        # Ambiguous base (N)
        elif base in 'Nn':
            base_quality = ord(qualities[q_pos]) - 33
            last_base = base
            last_base_quality = base_quality

            if base_quality >= min_bq:
                n_count += 1
            else:
                mm_count_lq += 1

            b_pos += 1
            q_pos += 1

        # End of read sign
        elif base == '$':
            b_pos += 1

        # Start of read sign, followed by MapQ ascii value
        elif base == '^':
            b_pos += 2

        # Deletion
        elif base == '-':

            # The . or , character before an indel is the indel anchor. Its quality represents the indels quality. It is not representing a ref call.
            if last_base_quality >= min_bq:
                ref_count = ref_count - 1
                if last_base == ".":
                    ref_count_fwd = ref_count_fwd - 1
                elif last_base == ",":
                    ref_count_rev = ref_count_rev - 1
            else:
                ref_count_lq = ref_count_lq - 1

            # One more deletion-supporting read found
            indel_counts['del'] += 1

            # Get deletion length and the deletion sequence
            b_pos += 1
            indel_len = bases[b_pos]

            while bases[b_pos + 1].isdigit():
                indel_len = indel_len + bases[b_pos + 1]
                b_pos += 1

            indel_len = int(indel_len)
            del_observed = bases[(int(b_pos) + 1):(int(b_pos) + int(indel_len) + 1)]
            strand = fwd_rev_strand(del_observed)
            del_observed = del_observed.upper()

            # Update count of previously observed deletion
            if del_observed in deletion:
                deletion[del_observed] += 1
                del_quals[del_observed] += last_base_quality

                if strand == "forward":
                    if del_observed in deletion_fwd:
                        deletion_fwd[del_observed] += 1
                    else:
                        deletion_fwd[del_observed] = 1

                elif strand == "reverse":
                    if del_observed in deletion_rev:
                        deletion_rev[del_observed] += 1
                    else:
                        deletion_rev[del_observed] = 1

            # Enter new deletion not observed before
            else:
                deletion[del_observed] = 1
                del_quals[del_observed] = last_base_quality

                if strand == "forward":
                    deletion_fwd[del_observed] = 1
                    deletion_rev[del_observed] = 0

                elif strand == "reverse":
                    deletion_rev[del_observed] = 1
                    deletion_fwd[del_observed] = 0

            # Jump over deletion sequence and continue to next base
            b_pos += int(indel_len)
            b_pos += 1

        # Insertion
        elif base == '+':

            # The . or , character before an indel is the indel anchor. Its quality represents the indels quality. It is not representing a ref call.
            if last_base == "." and last_base_quality >= min_bq:
                ref_count_fwd = ref_count_fwd - 1
                ref_count = ref_count - 1
            elif last_base == "," and last_base_quality >= min_bq:
                ref_count_rev = ref_count_rev - 1
                ref_count = ref_count - 1

            # One more insertion-supporting read found
            indel_counts['ins'] += 1

            # Get insertion length and insertion sequence
            b_pos += 1
            indel_len = bases[b_pos]

            while bases[b_pos + 1].isdigit():
                indel_len = indel_len + bases[b_pos + 1]
                b_pos += 1

            indel_len = int(indel_len)
            ins_observed = bases[(int(b_pos) + 1):(int(b_pos) + int(indel_len) + 1)]
            strand = fwd_rev_strand(ins_observed)
            ins_observed = ins_observed.upper()

            # Update count of previously observed insertion
            if ins_observed in insertion:
                insertion[ins_observed] += 1
                ins_quals[ins_observed] += last_base_quality

                if strand == "forward":
                    if ins_observed in insertion_fwd:
                        insertion_fwd[ins_observed] += 1
                    else:
                        insertion_fwd[ins_observed] = 1

                elif strand == "reverse":
                    if ins_observed in insertion_rev:
                        insertion_rev[ins_observed] += 1
                    else:
                        insertion_rev[ins_observed] = 1

            # Enter new insertion not observed before
            else:
                insertion[ins_observed] = 1
                ins_quals[ins_observed] = last_base_quality

                if strand == "forward":
                    insertion_fwd[ins_observed] = 1
                    insertion_rev[ins_observed] = 0
                elif strand == "reverse":
                    insertion_rev[ins_observed] = 1
                    insertion_fwd[ins_observed] = 0

            # Jump over insertion sequence and continue to next base
            b_pos += int(indel_len)
            b_pos += 1

        # Deletion placeholder. Deleted bases also have a quality character, which needs to be skipped.
        elif base == '*':
            b_pos += 1
            q_pos += 1

        # Error: unknown pileup character encountered
        else:
            print("ERROR:", locus, pos, b_pos, base)
            exit(1)

    # Compute depth and high-quality depth
    dp, dp_hq = coverage(qualities, min_bq)

    # Correct high quality depth: do not count Ns and indel-anchors
    dp_hq = dp_hq - n_count + indel_counts['ins'] + indel_counts['del']

    # Add reference counts to alt_counts dictionary
    alt_counts[ref_base.upper()] += ref_count_fwd
    alt_counts[ref_base.lower()] += ref_count_rev

    # Sum up counts from fwd and rev strand
    alt_count_total = {'A': alt_counts['A'] + alt_counts['a'], 'C': alt_counts['C'] + alt_counts['c'],
                       'G': alt_counts['G'] + alt_counts['g'], 'T': alt_counts['T'] + alt_counts['t']}

    # Sum up base qualities from fwd and rev strand
    alt_qual_total = {'A': alt_quals['A'] + alt_quals['a'], 'C': alt_quals['C'] + alt_quals['c'],
                      'G': alt_quals['G'] + alt_quals['g'], 'T': alt_quals['T'] + alt_quals['t']}

    # A counts
    change = str(ref_base.upper()) + '>A'
    a = [change, str(alt_count_total['A']), str(alt_counts['A'.upper()]), str(alt_counts['A'.lower()]),
         str(alt_qual_total['A']), '0']
    a = '\t'.join(a)

    # C counts
    change = str(ref_base.upper()) + '>C'
    c = [change, str(alt_count_total['C']), str(alt_counts['C'.upper()]), str(alt_counts['C'.lower()]),
         str(alt_qual_total['C']), '0']
    c = '\t'.join(c)

    # T counts
    change = str(ref_base.upper()) + '>T'
    t = [change, str(alt_count_total['T']), str(alt_counts['T'.upper()]), str(alt_counts['T'.lower()]),
         str(alt_qual_total['T']), '0']
    t = '\t'.join(t)

    # G counts
    change = str(ref_base.upper()) + '>G'
    g = [change, str(alt_count_total['G']), str(alt_counts['G'.upper()]), str(alt_counts['G'.lower()]),
         str(alt_qual_total['G']), '0']
    g = '\t'.join(g)

    # Insertion counts
    if indel_counts['ins'] > 0:
        ins_call = key_with_max_value(insertion)
        ins_count = insertion[ins_call]
        ins_count_fwd = insertion_fwd[ins_call]
        ins_count_rev = insertion_rev[ins_call]
        ref_base = ref_base.upper()
        ins_call = ref_base.upper() + ins_call.upper()
        diff = indel_counts['ins'] - ins_count
    else:
        ins_count = 0
        ins_count_fwd = 0
        ins_count_rev = 0
        diff = 0
        ref_base = ref_base.upper()
        ins_call = 'ins'

    change = ref_base + '>' + ins_call
    insertion = [change, str(ins_count), str(ins_count_fwd), str(ins_count_rev), '.', str(diff)]
    insertion = '\t'.join(insertion)

    # Deletion counts
    if indel_counts['del'] > 0:
        del_call = key_with_max_value(deletion)
        del_count = deletion[del_call]
        del_count_fwd = deletion_fwd[del_call]
        del_count_rev = deletion_rev[del_call]
        diff = indel_counts['del'] - del_count
        ref_base = ref_base.upper() + del_call.upper()
        del_call = ref_base.upper()
    else:
        del_count = 0
        del_count_fwd = 0
        del_count_rev = 0
        diff = 0
        ref_base = ref_base.upper()
        del_call = 'del'

    change = ref_base + '>' + del_call
    deletion = [change, str(del_count), str(del_count_fwd), str(del_count_rev), '.', str(diff)]
    deletion = '\t'.join(deletion)

    # Final result
    line = [str(locus), str(pos), str(ref_base), str(cov), str(dp_hq), str(ref_count_fwd), str(ref_count_rev),
            a, c, t, g, insertion, deletion, str(mm_count_lq), str(n_count)]
    line = '\t'.join(line)

    return line


# Run pileup_parser
def main():
    # Read parameters
    parser = argparse.ArgumentParser(description='Parse pileup file and prepare input tsv file for variant calling.')
    parser.add_argument('-p', '--pileup', type=str, help='Pileup file.', required=True)
    parser.add_argument('-o', '--outfile', type=str, required=True, dest='outfile',
                        help='Output file for predicted SNPs and indels.')
    parser.add_argument('-q', '--minBQ', required=False, dest='minBQ', type=int, default=10,
                        help='Minimum base quality to consider nucleotide for variant analysis. Default = 10')
    parser.add_argument('-m', '--minDP', required=False, dest='minDP', type=int, default=1,
                        help='Minimum total depth. Default = 1')

    try:
        args = parser.parse_args()
    except IOError as io:
        print(io)
        sys.exit('Error reading parameters.')

    input_pileup_file = args.pileup
    output_tsv_file = args.outfile
    min_bq = args.minBQ
    min_dp = args.minDP

    # Open output file to write the aggregated counts per genomic position
    output_tsv_file = open(output_tsv_file, 'w')

    # Header of the output tsv file
    header = (
        "CHROM", "POS", "REF", "DP", "DP_HQ", "REFf", "REFr", "Ag", "A", "Af", "Ar", "Aq", "Ao", "Cg", "C", "Cf", "Cr",
        "Cq", "Co", "Tg", "T", "Tf", "Tr", "Tq", "To", "Gg", "G", "Gf", "Gr", "Gq", "Go", "INSg", "INS", "INSf", "INSr",
        "INSq", "INSo", "DELg", "DEL", "DELf", "DELr", "DELq", "DELo", "Lq_alt_count", "N_count")
    header = ("\t".join(header)) + "\n"
    output_tsv_file.write(header)

    # Parse pileup file and send each row to the function pileup_info(line)
    with open(input_pileup_file) as f1:
        for pileup_line in f1:
            split_pileup_line = pileup_line.split("\t")

            # If pileup line has 6 columns and coverage is greater or equal to 5
            if len(split_pileup_line) == 6 and int(split_pileup_line[3]) >= min_dp:
                info = pileup_info(pileup_line, min_bq)
                output_tsv_file.write(info + "\n")

    # Finish
    output_tsv_file.close()
    exit(0)


main()
