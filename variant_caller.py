# coding=utf-8

# Import Libraries
import argparse
import time
import numpy
import pybedtools
import scipy.stats
from pathlib import Path

numpy.seterr(divide='ignore')


# Check sequence complexity (detects sequences with homopolymers)
def longest_run(s):
    if len(s) == 0:
        return 0

    runs = ''.join('*' if x == y else ' ' for x, y in zip(s, s[1:]))
    star_strings = runs.split()

    if len(star_strings) == 0:
        return 1

    return 1 + max(len(stars) for stars in star_strings)


# Get the percentage of the most frequent nucleotide in a sequence (identifies low complexity sequence)
def frequent_base(s):
    nucleotide_list = list(s)
    majority_base = max([nucleotide_list.count(base) for base in set(nucleotide_list)])
    allele_frequency = round(float(majority_base) / len(s), 2)
    return allele_frequency


# Obtain adjacent sequence
def up_down_sequence(chrom, start, up_down_length, infile, chrom_length=-1):
    # DEBUG: return arbitrary sequence. Done for speed check.
    sequence_list = ["ACGTG", "CGTAG"]
    return sequence_list


'''
    # As it works with 0-based coordinates, we must subtract 1 base to our start coordinate
    start = int(start)
    start = start - 1

    # Upstream
    end = start - up_down_length + 1

    # return 'N'-sequence if interval extends over chromosome bounds
    if int(end) < 1:
        seq_up = 'N' * up_down_length
    else:
        # Getting sequence downstream
        a = pybedtools.BedTool("\t".join([chrom, str(end), str(start)]), from_string=True)
        a = a.sequence(fi=infile)

        j = open(a.seqfn).read()
        j = j.rstrip('\n')

        # Sequence
        seq_up = j.split('\n')[1]

    # Downstream
    start = start + 1
    end = start + up_down_length

    # return 'N'-sequence if interval extends over chromosome bounds
    if chrom_length != -1 and end > chrom_length:
        seq_down = 'N' * up_down_length
    else:
        # Getting sequence upstream
        a = pybedtools.BedTool("\t".join([chrom, str(start), str(end)]), from_string=True)
        a = a.sequence(fi=infile)

        j = open(a.seqfn).read()
        j = j.rstrip('\n')
        seq_down = j.split('\n')[1]
    
    # Return upstream and downstream sequence    
    sequence_list = [seq_up, seq_down]
    return sequence_list
'''


# Filter variants using various filter types
def variant_filter(variant, mm, depth, depth_hq, dist, ref_t, seq_upstream, seq_downstream):
    gt, sig, alt_count, ab, alt_count_p, alt_count_fdr, strand_counts, strand_bias, alt_count_other, out_adj = variant.split(":")

    # Analysis of flanking sequences: check homopolymer length and sequence complexity
    seq_upstream_l = longest_run(seq_upstream)
    seq_upstream_freq = frequent_base(seq_upstream)

    seq_downstream_l = longest_run(seq_downstream)
    seq_downstream_freq = frequent_base(seq_downstream)

    # Reference and alternative base for this variant
    ref_base, alt_base = gt.split(">")

    # Get counts per strand
    alt_fwd, alt_rev, ref_fwd, ref_rev = strand_counts.split("-")

    # Set filter based on various thresholds
    if float(alt_count) > 0 and mm > 0:
        filter_criteria = []

        if float(alt_count_fdr) > 0.1:
            filter_criteria.append("Error")

        if (seq_upstream_l >= 3 or seq_upstream_freq >= 0.8) and len(gt) > 3:
            filter_criteria.append("LC_Upstream")

        if (seq_downstream_l >= 3 or seq_downstream_freq >= 0.8) and len(gt) > 3:
            filter_criteria.append("LC_Downstream")

        if float(ab) < min_AF:
            filter_criteria.append("Low_AF")

        if (min(int(alt_fwd), int(alt_rev)) == 0 and
                min(int(ref_fwd), int(ref_rev)) != 0 and
                strand_counts == 1 and int(alt_count) > 3):
            filter_criteria.append("Strand_imbalanced")

        if int(depth_hq) < min_COV:
            filter_criteria.append("Low_Cov")

        if int(alt_count) < int(min_AC):
            filter_criteria.append("Low_AC")

        if dist != "Inf" and int(dist) < min_DIST:
            filter_criteria.append("Clustered_Variant")

        if float(depth_hq) / float(depth) < 0.75:
            filter_criteria.append("Low_qual_pos")

        if float(alt_count_other) / (float(mm)) > 0.3 or float(out_adj) < 0.1:
            filter_criteria.append("Variant_contamination")

        if float(strand_bias) < 0.01 and strand_counts == 1:
            filter_criteria.append("Fisher_Strand")

        if len(filter_criteria) == 0:
            concatenated_filter_string = "PASS"
        else:
            concatenated_filter_string = ';'.join(filter_criteria)
    else:
        alt_base = '.'
        concatenated_filter_string = '.'

    tsv = [alt_count, depth_hq, ab, ref_t, alt_count_p, alt_count_fdr, strand_counts, strand_bias, alt_count_other, out_adj,
           str(seq_upstream_freq), str(seq_downstream_freq)]

    variant_call = [ref_base, alt_base, tsv, concatenated_filter_string]

    return variant_call


# Print VCF header and column names
def print_vcf_header():

    # Define VCF header
    date = time.strftime("%d/%m/%Y")  # dd/mm/yyyy format
    vcf_format = "##fileformat=VCFv4.1"
    date_line = "##fileDate=%s" % date
    source = "##source=CRG_UKT_somatic_variant_calling"
    reference = "##reference=%s" % reference_file
    concepts = """##INFO=<ID=Variant_Dist,Number=1,Type=Integer,Description="Distance to the closest short variant">
    ##INFO=<ID=Upstream,Number=.,Type=String,Description="Upstream sequence (5 nucleotides)">
    ##INFO=<ID=Downstream,Number=.,Type=String,Description="Downstream sequence (5 nucleotides)">
    ##INFO=<ID=PValue,Number=.,Type=String,Description="Uncorrected p-value">
    ##FILTER=<ID=PASS,Description="Passed filter">
    ##FILTER=<ID=Low_COV,Description="Low coverage">
    ##FILTER=<ID=Strand_imbalanced,Description="All alternative reads found in only one strand">
    ##FILTER=<ID=Low_AC,Description="Less than defined minimum of alternative counts">
    ##FILTER=<ID=Clustered_Variant,Description="Clustered variants">
    ##FILTER=<ID=LC_Upstream,Description="Low complexity region (5bps) upstream. ≥ 80% of bases show the same nucleotide or tandem of ≥ 3 equal nucleotides in a row">
    ##FILTER=<ID=LC_Downstream,Description="Low complexity region (5bps) downstream. ≥ 80% of bases show the same nucleotide  or tandem of ≥ 3 equal nucleotides in a row">
    ##FILTER=<ID=Error,Description="Alternative counts inside the expected error rate distribution">
    ##FILTER=<ID=Fisher_Strand,Description="Strand bias based on fisher test">
    ##FILTER=<ID=Low_qual_pos,Description="Position enriched with too many low quality bases">
    ##FILTER=<ID=Variant_contamination,Description="Reads supporting other alleles outsite of the error rate distribution">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele read counts"
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance">
    ##FORMAT=<ID=Strand,Number=2,Type=String,Description="Alleles in strands: Alternative forward, Alternative reverse, Reference forward, Reference reverse">
    ##FORMAT=<ID=FS,Number=1,Type=Float,Description="Fisher strand test (q-value)">
    ##FORMAT=<ID=VCB,Number=1,Type=Integer,Description="Variant Count bias, number of other different alternative alleles found">
    ##FORMAT=<ID=Perror,Number=1,Type=Float,Description="Q-value to belong to the error rate distribution (beta binomial distribution) - After FDR correction">"""

    # Write VCF header and column names
    vcf_column_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_id]
    vcf_column_names_line = "#" + '\t'.join(vcf_column_names)
    OUT_vcf.write(vcf_format + '\n')
    OUT_vcf.write(date_line + '\n')
    OUT_vcf.write(source + '\n')
    OUT_vcf.write(reference + '\n')
    OUT_vcf.write(concepts + '\n')
    OUT_vcf.write(vcf_column_names_line + '\n')


# Write column names to tsv outfile
def print_tsv_header():
    tsv_header = ['CHROM', 'POS', 'REF', 'ALT', 'Upstream_5', 'Downstream_5', 'DP_HQ', 'REFt', 'ALT_COUNT', 'AB', 'P_VAL',
                  'P_VAL_adj', 'STRAND', 'FISHER', 'ALT_COUNT_o', 'P_VALo_adj', 'LC_Upstream', 'LC_Downstream', 'FILTER']
    tsv_header = '\t'.join(tsv_header)
    OUT_tsv.write(tsv_header + '\n')


# Parse aggregated pileup statistics (tsv) file with all regions of interest and call variants
def parse_pileup_statistics():

    # Parse aggregated pileup statistics (tsv) file with all regions of interest and call variants
    with open(input_file) as f1:
        for line in f1:
            line = line.rstrip('\n')

            # Ignore header lines
            if line.startswith('##'):
                continue

            # Read column names and content of input file
            elif line.startswith('CHROM'):
                header_file = line
                form = header_file.split('\t')

                chrom_index = [i for i, x in enumerate(form) if x == "CHROM"][0]
                pos_index = [i for i, x in enumerate(form) if x == "POS"][0]
                dp_index = [i for i, x in enumerate(form) if x == "DP"][0]
                dp_hq_index = [i for i, x in enumerate(form) if x == "DP_HQ"][0]
                ref_fwd_index = [i for i, x in enumerate(form) if x == "REFf"][0]
                ref_rev_index = [i for i, x in enumerate(form) if x == "REFr"][0]
                dist_index = [i for i, x in enumerate(form) if x == "DIST"][0]
                mm_index = [i for i, x in enumerate(form) if x == "MM"][0]
                all_index = [i for i, x in enumerate(form) if x == "CALL"][0]

            # Parse rows with actual loci data
            else:
                info_field = line.split("\t")

                # Assign values to common VCF columns
                chrom = info_field[chrom_index]
                pos = info_field[pos_index]
                snv_id = '.'
                ref_fwd = int(info_field[ref_fwd_index])
                ref_rev = int(info_field[ref_rev_index])
                ref_total = ref_fwd + ref_rev
                mm = int(info_field[mm_index])
                dp = int(info_field[dp_index])
                dp_hq = int(info_field[dp_hq_index])
                qual = '.'
                dist = info_field[dist_index]

                # Determine length of the current chromosome
                chr_length = -1
                if chr_lengths != {}:
                    if chrom in chr_lengths.keys():
                        chr_length = chr_lengths[chrom]

                # Getting 5 bases up and downstream of focal locus
                seq_up, seq_down = up_down_sequence(chrom, pos, 6, args.reference, chr_length)

                # Read count info
                call = info_field[all_index].split("|")

                # print CHROM, POS
                ref = []
                alt = []
                tsv = []
                filter_string = []
                for variant in call:
                    variant_ref, variant_alt, variant_tsv, variant_filter_string = variant_filter(variant, mm, dp, dp_hq, dist, ref_total, seq_up, seq_down)

                    # Append variants
                    ref.append(variant_ref)
                    alt.append(variant_alt)
                    tsv.append(variant_tsv)
                    filter_string.append(variant_filter_string)

                if len(ref) > 1:

                    # Getting ref. For both SNP, Deletion and Insertion, the longer reference will represent the original reference
                    reference_seq = max(ref, key=len)

                    filter_field = []
                    alt_base = []
                    alt_count_all = []
                    sample = []

                    n = len(ref)

                    p_val_list = []
                    for count in range(0, n):
                        ref_index = ref[count]
                        if len(ref_index) != len(reference_seq):
                            # Sample info
                            alt_count, dp_hq, ab_score, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, seq_upstream_freq, seq_downstream_freq = tsv[count]

                            # Getting normalized alternative allele
                            alt_index = alt[count]
                            alt_temp = list(reference_seq)
                            alt_temp[0] = alt_index
                            alt_index = ''.join(alt_temp)
                            alt_base.append(alt_index)
                            alt_count_all.append(alt_count)

                            if int(alt_count) > 0:
                                genotype = alt_index
                            else:
                                genotype = reference_seq

                            filter_field.append(filter_string[count])

                            sample_index = [genotype, str(alt_count), str(dp_hq), str(ab_score), base_strand, str(fisher), str(alt_count_o), str(out_adj), str(alt_count_pad_j)]
                            sample.append(':'.join(sample_index))

                            p_val_list.append(float(alt_count_p))

                        else:

                            # Sample info
                            alt_count, dp_hq, ab_score, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, seq_upstream_freq, seq_downstream_freq = tsv[count]

                            alt_index = alt[count]
                            alt_base.append(alt_index)
                            alt_count_all.append(alt_count)

                            if int(alt_count) > 0:
                                genotype = alt_index
                            else:
                                genotype = reference_seq

                            filter_field.append(filter_string[count])

                            sample_index = [genotype, str(alt_count), str(dp_hq), str(ab_score), base_strand, str(fisher), str(alt_count_o), str(out_adj), str(alt_count_pad_j)]
                            sample.append(':'.join(sample_index))

                            p_val_list.append(float(alt_count_p))

                    # Get number of alternative alleles
                    if len(alt_base) != len(filter_field) or len(alt_base) != len(sample) or len(alt_base) != len(alt_count_all):
                        raise ValueError("Number of alt alleles differs")

                    # Combine p-values of this site
                    p_merge = scipy.stats.combine_pvalues(p_val_list)[1]
                    MRD.append(float(p_merge))

                    # Split multi-allelic variants in single lines
                    for allele_idx in range(len(alt_base)):
                        # Common vcf columns
                        vcf_fields = [chrom, pos, snv_id, reference_seq, alt_base[allele_idx], qual]

                        # Format column
                        format_field = ["GT", "Alt_Count", "DP", "AB", "Strand", "FS", "VCB", "Pvcb", "Perror"]

                        # Info column
                        info_field = ['Variant_dist=' + str(dist), 'Upstream=' + str(seq_up), 'Downstream=' + str(seq_down), 'PValue=' + str(alt_count_p)]

                        # VCF variant line
                        vcf_line = ['\t'.join(vcf_fields), filter_field[allele_idx], ';'.join(info_field), ':'.join(format_field),
                                    sample[allele_idx]]
                        vcf_line = '\t'.join(vcf_line)

                        # TSV line
                        tsv_line = [str(chrom), str(pos), str(reference_seq), str(alt_base[allele_idx]), seq_up, seq_down, str(dp_hq), str(ref_total),
                                    str(alt_count_all[allele_idx]), str(ab_score), alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o,
                                    out_adj, seq_upstream_freq, seq_downstream_freq, filter_field[allele_idx]]
                        tsv_line = '\t'.join(tsv_line)

                        if int(alt_count) > 0:
                            OUT_vcf.write(vcf_line + '\n')
                            OUT_tsv.write(tsv_line + '\n')

                else:

                    # Common vcf columns
                    reference_seq = ref[0]
                    alt_base = alt[0]
                    vcf_fields = [chrom, pos, snv_id, reference_seq, alt_base, qual]

                    # Filter column
                    filter_field = filter_string[0]

                    # Sample info
                    alt_count, dp_hq, ab_score, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, seq_upstream_freq, seq_downstream_freq = tsv[0]

                    # Info column
                    info_field = ['Variant_dist=' + str(dist), 'Upstream=' + str(seq_up), 'Downstream=' + str(seq_down), 'PValue=' + str(alt_count_p)]

                    # Format column
                    format_field = ["GT", "Alt_Count", "DP", "AB", "Strand", "FS", "VCB", "Pvcb", "Perror"]

                    # Append p-val for MRD calculation
                    MRD.append(float(alt_count_p))

                    # Get genotype
                    if int(alt_count) > 0:
                        genotype = alt_base
                    else:
                        genotype = reference_seq

                    sample = [genotype, str(alt_count), str(dp_hq), str(ab_score), base_strand, str(fisher), str(alt_count_o), str(out_adj),
                              str(alt_count_pad_j)]

                    # VCF variant line
                    vcf_line = ['\t'.join(vcf_fields), filter_field, ';'.join(info_field), ':'.join(format_field), ':'.join(sample)]
                    vcf_line = '\t'.join(vcf_line)

                    # TSV line
                    tsv_line = [str(chrom), str(pos), str(ref[0]), str(alt[0]), seq_up, seq_down, str(dp_hq), str(ref_total),
                                str(alt_count), str(ab_score), alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj,
                                seq_upstream_freq, seq_downstream_freq, filter_field]
                    tsv_line = '\t'.join(tsv_line)
                    if int(alt_count) > 0:
                        OUT_vcf.write(vcf_line + '\n')
                        OUT_tsv.write(tsv_line + '\n')


# Compute and write Minimal Residual Disease probability
def compute_mrd():
    out3 = str(Path(vcf_outfile_name).with_suffix('.mrd'))
    out_mrd = open(out3, 'w')

    # Header
    mrd_header = ['#MRD_log10', 'MRD_pval']
    mrd_header = '\t'.join(mrd_header)
    out_mrd.write(mrd_header + '\n')

    # Getting MRD values
    numpy.seterr(divide='ignore')
    val_mrd_temp = scipy.stats.combine_pvalues(MRD)[1]
    val_mrd_temp2 = max(1e-20, val_mrd_temp)
    val_mrd = str(format(val_mrd_temp2, "5.2e"))
    val_mrd_log = str(round(numpy.log10(float(val_mrd)) * -1, 4))

    # Printing line
    line = [val_mrd_log, val_mrd]
    line = '\t'.join(line)
    out_mrd.write(line + '\n')
    out_mrd.close()


# Main method

# Parse arguments
parser = argparse.ArgumentParser(description='Getting barcodes in fastq file and labels')
parser.add_argument('-i', '--infile', type=str, help='Tsv table', required=True)
parser.add_argument('-tID', '--tumorid', type=str, default='Tumor', help='Tumor sample id', required=False)
parser.add_argument('-ref', '--reference', type=str, help='Reference fasta file which table was build', required=True)
parser.add_argument('-o', '--outfile', type=str, help='Vcf output file', required=True)
parser.add_argument('-cov', '--min_COV', type=int, default=10, help='Minimum Coverage', required=False)
parser.add_argument('-ac', '--min_AC', type=int, default=3, help='Minimum reads supporting alternative allele', required=False)
parser.add_argument('-variant_dist', '--min_DIST', type=int, default=20, help='Minimum distance allowed between variants (to avoid clustered errors)', required=False)
parser.add_argument('-str', '--strand', type=int, choices=[0, 1], default=1, help='Strand bias test (Fisher test). 0 for turn it off', required=False)
parser.add_argument('-af', '--min_AF', type=float, default=0, help='Minimum allele frequency allowed', required=False)
parser.add_argument('-mrd', '--mrd', type=int, choices=[0, 1], default=1, help='Print Minimal Residual Disease [Default = 1]', required=False)
parser.add_argument('-tmpdir', '--tmpdir', default=None, help='Folder for temp files', required=False)


# Store arguments
args = parser.parse_args()
input_file = args.infile
sample_id = args.tumorid
reference_file = args.reference
min_COV = args.min_COV
min_AC = args.min_AC
min_AF = args.min_AF
min_DIST = args.min_DIST
strand = args.strand
MRD = []


# Check temp folder
if args.tmpdir is not None:
    pybedtools.set_tempdir(args.tmpdir)


# Open output files
vcf_outfile_name = args.outfile
tsv_outfile_name = str(Path(vcf_outfile_name).with_suffix('.tsv'))
OUT_vcf = open(vcf_outfile_name, 'w')
OUT_tsv = open(tsv_outfile_name, 'w')


# Load reference Fasta index to determine chr lengths
chr_lengths = {}
if Path(args.reference + ".fai").exists():
    with open(args.reference + ".fai", "rt") as fasta_index_file:
        for fasta_line in fasta_index_file.readlines():
            split_line = fasta_line.split('\t')
            chr_lengths[split_line[0].strip()] = int(split_line[1].strip())
else:
    print("WARNING: No FASTA index file found! Can't determine chromosome lengths.")


# Print output headers and column names
print_vcf_header()
print_tsv_header()


# Parse the input file:
parse_pileup_statistics()


# Compute Minimal Residual Disease
if args.mrd == 1:
    compute_mrd()


# Close output files
OUT_vcf.close()
OUT_tsv.close()
exit(0)
