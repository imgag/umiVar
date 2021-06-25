# coding=utf-8

# Import Libraries
import argparse
import os
import time
import numpy
import scipy.stats
from pathlib import Path
import pysam

numpy.seterr(divide='ignore')


# Returns longest homopolymer
def longest_homopolymer(sequence):
    if len(sequence) == 0:
        return 0

    runs = ''.join('*' if x == y else ' ' for x, y in zip(sequence, sequence[1:]))
    star_strings = runs.split()
    if len(star_strings) == 0:
        return 1

    return 1 + max(len(stars) for stars in star_strings)


# Returns most frequent base in a string
def frequent_base(sequence):
    # Get the most common base and its percentage of length of the sequence
    nucleotide_list = list(sequence)
    most_common_base = max([nucleotide_list.count(base) for base in set(nucleotide_list)])
    most_common_base_percentage = round(float(most_common_base) / len(sequence), 2)

    return most_common_base_percentage


# Homopolymer filter
def homopolymer_filter(sequence):
    base_type_zero_count = 0  # number of base types (A, C, G, T) not represented in the sequence
    low_complexity_flag = 0  # low complexity = (less than 3 base types represented or homopolymers of length >= 4)
    homopolymer_flag = 0  # homopolymer of length >= 5

    if sequence != '.':

        # Get the longest k-mer
        max_homopolymer_length = longest_homopolymer(sequence)

        # Get the frequency of the homopolymer base
        max_base_freq = frequent_base(sequence)

        # Count base types
        nucleotide_list = list(sequence)
        for base in ['A', 'C', 'G', 'T']:
            if nucleotide_list.count(base) == 0:
                base_type_zero_count += 1

        # Set homopolymer flag
        if max_homopolymer_length >= 5 or max_base_freq >= 0.8:
            homopolymer_flag = 1

        # Set low complexity flag
        if base_type_zero_count > 1 or max_homopolymer_length >= 4:
            low_complexity_flag = 1

    return homopolymer_flag, low_complexity_flag


# Filter variants using various filter types
def variant_filter(variant, mm, depth, depth_hq, dist, ref_t, sequence_context):
    # Declare variables
    gt, sig, alt_count, af, alt_p, alt_fdr, strand_counts, strand_bias, other_alt, other_fdr = variant.split(":")

    # Reference and alternative base for this variant
    ref_base, alt_base = gt.split(">")

    # Get counts per strand
    alt_fwd, alt_rev, ref_fwd, ref_rev = strand_counts.split("-")

    # Homopolymer filter on upstream / downstream region
    homopolymer, low_complexity = homopolymer_filter(sequence_context)

    # Set filter based on various thresholds
    if float(alt_count) > 0 and mm > 0:
        filter_criteria = []

        # Filter: p-value not significant
        if float(alt_p) > 0.1:
            filter_criteria.append("pvalue")

        # Filter: FDR not significant
        if float(alt_fdr) > 0.1:
            filter_criteria.append("FDR")

        # Filter: homopolymer
        if homopolymer == 1:
            filter_criteria.append("Homopolymer")

        # Filter: low complexity
        if low_complexity == 1:
            filter_criteria.append("Low_Complexity")

        # Filter: minimum allele frequency
        if float(af) < min_AF:
            filter_criteria.append("Low_AF")

        # Strand imbalanced (very low alt counts on one strand)
        if int(alt_fwd) <= 1 or int(alt_rev) <= 1:
            filter_criteria.append("Strand_Imbalanced")

        # Filter: minimum depth
        if int(depth_hq) < min_COV:
            filter_criteria.append("Low_Cov")

        # Minimum alternative allele count filter
        if int(alt_count) < int(min_AC):
            filter_criteria.append("Low_AC")

        # Filter: clustered variants
        if dist != "Inf" and int(dist) < min_DIST:
            filter_criteria.append("Clustered_Variant")

        # Filter: High fraction of low quality bases
        if float(depth_hq) / float(depth) < 0.75:
            filter_criteria.append("Low_Qual_Pos")

        # Filter: contamination = high number of third or fourth allele
        if float(other_alt) / (float(mm)) > 0.3 or float(other_fdr) < 0.1:
            filter_criteria.append("Variant_Contamination")

        # Filter: Fisher strand bias
        if float(strand_bias) < 0.01 and strand_counts == 1:
            filter_criteria.append("Fisher_Strand")

        # Final filter judgement (PASS if no filter criteria apply)
        if len(filter_criteria) == 0:
            concatenated_filter_string = "PASS"
        else:
            concatenated_filter_string = ';'.join(filter_criteria)
    else:
        alt_base = '.'
        concatenated_filter_string = '.'

    # Set output columns for tsv format
    tsv = [alt_count, depth_hq, af, ref_t, alt_p, alt_fdr, strand_counts, strand_bias, other_alt, other_fdr, homopolymer]

    # Compile and return list with results
    variant_call = [ref_base, alt_base, tsv, concatenated_filter_string]
    return variant_call


# Print variants of various quality levels to separate files
def print_variants(vcf_line, tsv_line):

    # Extract fields of tsv entry
    tsv_entries = tsv_line.split("\t")

    # Check if printing is required (at least one alternative read)
    ref_base = tsv_entries[2]
    alt_base = tsv_entries[3]
    alt_count = tsv_entries[6]
    multi_umi_alt_depth = tsv_entries[9]
    filter_string = tsv_entries[19]

    # Check if indel
    is_indel = False
    if len(ref_base) > 1 or len(alt_base) > 1:
        is_indel = True

    # Check if high quality variant
    high_quality = True

    # General quality parameters
    if "Homopolymer" in filter_string or \
            "Low_AF" in filter_string or \
            "Low_AC" in filter_string or \
            "pvalue" in filter_string or \
            "Strand_Imbalanced" in filter_string:
        high_quality = False

    # Indel-specific quality check
    if is_indel:
        if "Clustered_Variant" in filter_string or "Low_Complexity" in filter_string:
            high_quality = False

    # Check if printing to verbose output file is required (at least one alternative multi-UMI read)
    if int(alt_count) > 0 and int(multi_umi_alt_depth) > 0:
        OUT_vcf.write(vcf_line + '\n')
        OUT_tsv.write(tsv_line + '\n')

        # Write high quality calls
        if high_quality:
            OUT_vcf_hq.write(vcf_line + '\n')
            OUT_tsv_hq.write(tsv_line + '\n')


# Parse aggregated pileup statistics (tsv) file with all regions of interest and call variants ------------------------
def parse_pileup_statistics():

    # Open fasta file with pysam
    in_fasta = pysam.FastaFile(reference_file)

    # Parse tsv file with counts of singleton UMIs (dedup_DP1.tsv) to get singleton counts ----------------------------
    singleton_dp = dict()
    with open(singleton_umi_file) as s1:
        for singleton_line in s1:

            # Ignore header lines
            if singleton_line.startswith('CHROM'):
                continue

            # Extract columns
            singleton_line = singleton_line.rstrip('\n')
            singleton_column = singleton_line.split("\t")

            # Generate locus-ID for dictionary
            locus_id = singleton_column[0] + "_" + singleton_column[1]

            # Store singleton-depth info in nested dictionary
            if locus_id not in singleton_dp:
                singleton_dp[locus_id] = dict()

            singleton_dp[locus_id]['A'] = singleton_column[8]
            singleton_dp[locus_id]['C'] = singleton_column[14]
            singleton_dp[locus_id]['T'] = singleton_column[20]
            singleton_dp[locus_id]['G'] = singleton_column[26]
            singleton_dp[locus_id]['I'] = singleton_column[32]
            singleton_dp[locus_id]['D'] = singleton_column[38]
            singleton_dp[locus_id]['.'] = 0

    # Parse input file stats.txt --------------------------------------------------------------------------------------
    with open(input_file) as f1:
        for line in f1:

            monitoring_flag = 0     # check if monitoring or ID position has been successfuly printed

            line = line.rstrip('\n')

            # Ignore header lines
            if line.startswith('##') or line.startswith('CHROM'):
                continue

            # Parse rows with actual data from stats.txt: CHROM, POS, REF, DP, DP_HQ, REFf, REFr, MM, CALL, DIST
            stats_fields = line.split("\t")

            # Assign values to common VCF columns
            chrom = stats_fields[0]
            pos = stats_fields[1]
            ref_nucleotide = stats_fields[2]
            dp = int(stats_fields[3])
            dp_hq = int(stats_fields[4])
            ref_fwd = int(stats_fields[5])
            ref_rev = int(stats_fields[6])
            mm = int(stats_fields[7])
            dist = stats_fields[9]
            ref_total = ref_fwd + ref_rev
            snv_id = '.'
            qual = '.'

            # Extract all possible alternative alleles per position:
            # Could include up to 3 nucleotide changes and multiple indel lengths -------------------------------------
            variant_candidates = stats_fields[8].split("|")

            # Get sequence context (vicinity) of a variant for homopolymer check (5 bases up- and down-stream)
            # Get fewer bases when variant is at the start or end of the sequence
            num_bases = 5
            pos_start = max(int(pos) - (num_bases + 1), 1)
            pos_end = min(int(pos) + num_bases, in_fasta.get_reference_length(chrom))
            try:
                sequence_context = in_fasta.fetch(chrom, pos_start, pos_end)
                sequence_context = sequence_context.upper()
            except FileNotFoundError:
                sequence_context = '.'

            # Containers for information about variant positions ------------------------------------------------------
            ref = []
            alt = []
            tsv = []
            filter_string = []

            # Extract all possible variants at the position and apply variant filter method to each possible variant
            for variant in variant_candidates:
                variant_ref, variant_alt, variant_tsv, variant_filter_string = variant_filter(variant, mm, dp, dp_hq, dist, ref_total, sequence_context)

                # Store variant in container
                ref.append(variant_ref)
                alt.append(variant_alt)
                tsv.append(variant_tsv)
                filter_string.append(variant_filter_string)

            # Go through all possible variants at the position --------------------------------------------------------
            for variant_iterator in range(0, len(alt)):

                # Sample info
                alt_count, dp_hq, allele_frequency, ref_total, alt_p, alt_fdr, base_strand, fisher, other_alt, other_fdr, homopolymer = tsv[variant_iterator]

                # Common vcf columns
                vcf_fields = [chrom, pos, snv_id, ref[variant_iterator], alt[variant_iterator], qual]

                # Filter column
                filter_field = filter_string[variant_iterator]

                # VCF Info column
                info_fields = ['Variant_dist=' + str(dist), 'Vicinity=' + str(sequence_context)]

                # VCF Format column
                format_field = ["GT", "DP", "AC", "AF", "M_REF", "M_AC", "M_AF", "Pval", "FDR", "Strand", "FS", "OAAC", "MAL", ]

                # Determine genotype:
                if float(allele_frequency) > 0.8:
                    genotype = "1/1"
                elif float(allele_frequency) > 0.3:
                    genotype = "0/1"
                else:
                    genotype = "./."

                # Extract count of singleton UMI fragments for alternative base at locus
                locus_id = str(chrom) + "_" + str(pos)
                alt_key = alt[variant_iterator]
                if len(ref[variant_iterator]) > 1:
                    alt_key = 'D'
                if len(alt_key) > 1:
                    alt_key = 'I'

                singleton_ref_depth = singleton_dp[locus_id][ref_nucleotide]
                singleton_alt_depth = singleton_dp[locus_id][alt_key]
                multi_umi_ref_depth = int(ref_total) - int(singleton_ref_depth)
                multi_umi_alt_depth = int(alt_count) - int(singleton_alt_depth)
                try:
                    multi_umi_af = multi_umi_alt_depth / (multi_umi_ref_depth + multi_umi_alt_depth)
                except ZeroDivisionError:
                    multi_umi_af = 0

                # VCF Sample column
                sample = [genotype, str(alt_count), str(dp_hq), str(allele_frequency),
                          str(multi_umi_ref_depth), str(multi_umi_alt_depth), str(round(multi_umi_af, 6)),
                          str(alt_p), str(alt_fdr), base_strand, str(fisher), str(other_alt), str(other_fdr)]

                # Compile VCF entry
                vcf_line = ['\t'.join(vcf_fields), filter_field, ';'.join(info_fields), ':'.join(format_field), ':'.join(sample)]
                vcf_line = '\t'.join(vcf_line)

                # Compile TSV entry
                tsv_line = [str(chrom), str(pos), str(ref[variant_iterator]), str(alt[variant_iterator]),
                            str(dp_hq), str(ref_total), str(alt_count), str(allele_frequency),
                            str(multi_umi_ref_depth), str(multi_umi_alt_depth), str(round(multi_umi_af, 6)),
                            alt_p, alt_fdr, base_strand, fisher, other_alt, other_fdr,
                            str(sequence_context), str(homopolymer), filter_field]
                tsv_line = '\t'.join(tsv_line)

                # Print variant candidates
                print_variants(vcf_line, tsv_line)

                # Print monitoring SNVs -------------------------------------------------------------------------------
                if bool(monitoring_variants):

                    # Check if position is in monitoring dictionary
                    if chrom in monitoring_variants.keys() and pos in monitoring_variants[chrom].keys():

                        # Setting monitoring flag to 1 indicates that this position must be printed, even if zero alts
                        if monitoring_flag == 0:
                            monitoring_flag = 1

                        # Check if base change matches monitoring or ID variant
                        if monitoring_variants[chrom][pos][1] == ref[variant_iterator] and \
                                monitoring_variants[chrom][pos][2] == alt[variant_iterator]:

                            # Setting monitoring flag to 1 indicates that position has been printed
                            monitoring_flag = 2

                            # Write either ID or monitoring variant
                            if monitoring_variants[chrom][pos][0] == 'ID':
                                TSV_id.write(tsv_line + '\n')
                                VCF_id.write(vcf_line + '\n')
                            elif monitoring_variants[chrom][pos][0] == 'M':
                                TSV_monitoring.write(tsv_line + '\n')
                                VCF_monitoring.write(vcf_line + '\n')

                                # Store p-val for MRD calculation
                                MRD_P.append(float(alt_p))
                                MRD_DP.append(int(ref_total) + int(alt_count))
                                MRD_ALT.append(int(alt_count))

            # Check if a monitoring position has not been printed
            if monitoring_flag == 1:

                # Compile VCF output
                vcf_fields = [chrom, pos, snv_id, str(monitoring_variants[chrom][pos][1]), str(monitoring_variants[chrom][pos][2]), qual]
                sample = ["0/0", "0", str(dp_hq), "0", str(multi_umi_ref_depth), "0", "0", "1.0", "1.0", "NA", "1.0", "0", "1.0"]
                vcf_line = ['\t'.join(vcf_fields), "NA", ';'.join(info_fields), ':'.join(format_field), ':'.join(sample)]
                vcf_line = '\t'.join(vcf_line)

                # Compile TSV output
                tsv_line = [str(chrom), str(pos), str(monitoring_variants[chrom][pos][1]), str(monitoring_variants[chrom][pos][2]),
                            str(dp_hq), str(ref_total), "0", "0", str(multi_umi_ref_depth), "0", "0", "1.0", "1.0",
                            "NA", "1.0", "0", "1.0", str(sequence_context), str(homopolymer), "NA"]
                tsv_line = '\t'.join(tsv_line)

                # Print monitoring positions that are homozygous reference
                if monitoring_variants[chrom][pos][0] == 'ID':
                    TSV_id.write(tsv_line + '\n')
                    VCF_id.write(vcf_line + '\n')
                elif monitoring_variants[chrom][pos][0] == 'M':
                    TSV_monitoring.write(tsv_line + '\n')
                    VCF_monitoring.write(vcf_line + '\n')
                    MRD_P.append(1.0)
                    MRD_DP.append(int(ref_total))
                    MRD_ALT.append(0)


# Compute and write Minimal Residual Disease probability
def compute_mrd():
    mrd_outfile = str(Path(vcf_outfile_name).with_suffix('.mrd'))
    mrd_fh = open(mrd_outfile, 'w')

    # Header
    mrd_fh.write('#MRD_log10' + '\t' + 'MRD_pval' + '\t' + 'SUM_DP' + '\t' + 'SUM_ALT' + '\t' + 'Mean_AF' + '\t' + 'Median_AF' + '\n')

    # Compute MRD pvalue and log10 transformation
    numpy.seterr(divide='ignore')
    val_mrd = scipy.stats.combine_pvalues(MRD_P)[1]
    val_mrd = max(1e-20, val_mrd)
    val_mrd = format(val_mrd, "5.2e")
    val_mrd_log = round(numpy.log10(float(val_mrd)) * -1, 4)

    # Compute aggregated allele fraction of monitoring mutations
    sum_dp = sum(MRD_DP)
    sum_alt = sum(MRD_ALT)
    print(str(sum_dp) + "\t" + str(sum_alt) + "\n")
    aggregated_af = float(sum_alt) / (float(sum_dp))

    # Compute allele fraction of medians
    median_af = float(numpy.median(MRD_ALT)) / float(numpy.median(MRD_DP))

    # Print results and close file
    mrd_fh.write(str(val_mrd_log) + '\t' + str(val_mrd) + '\t' + str(sum_dp) + "\t" + str(sum_alt) + "\t" + str(round(aggregated_af, 6)) + '\t' + str(round(median_af, 6)) + '\n')
    mrd_fh.close()


# Print VCF header and column names
def print_vcf_header():

    # Define VCF header
    date = time.strftime("%d/%m/%Y")  # dd/mm/yyyy format
    vcf_format = "##fileformat=VCFv4.1"
    date_line = "##fileDate=%s" % date
    source = "##source=umiVar2"
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
##FILTER=<ID=Variant_contamination,Description="Reads supporting other alleles outside of the error rate distribution">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele read counts">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternative allele fraction">
##FORMAT=<ID=M_REF,Number=1,Type=Integer,Description="Multi-UMI reference-like read counts">
##FORMAT=<ID=M_AC,Number=.,Type=Integer,Description="Multi-UMI alternative read counts">
##FORMAT=<ID=M_AF,Number=1,Type=Float,Description="Multi-UMI alternative allele fraction">
##FORMAT=<ID=Pval,Number=1,Type=Float,Description="P-value to belong to the error rate distribution (beta binomial distribution)">
##FORMAT=<ID=FDR,Number=1,Type=Float,Description="Q-value to belong to the error rate distribution (beta binomial distribution) after Benjamini-Hochberg correction (FDR)">
##FORMAT=<ID=Strand,Number=1,Type=String,Description="Alleles in strands: Alternative forward, Alternative reverse, Reference forward, Reference reverse">
##FORMAT=<ID=FS,Number=1,Type=Float,Description="Fisher strand bias test (q-value)">
##FORMAT=<ID=OAAC,Number=1,Type=Integer,Description="Other Alternative Allele Count = count of reads showing a different alternative allele">
##FORMAT=<ID=MAL,Number=1,Type=Float,Description="Multi-allele likelihood = Q-value (FDR) for other alternative alleles to belong to the error rate distribution (beta binomial distribution)">"""

    # Write VCF header and column names
    vcf_column_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_id]
    vcf_column_names_line = "#" + '\t'.join(vcf_column_names)

    # Compile header
    final_header = vcf_format + '\n' + date_line + '\n' + source + '\n' + reference + '\n' + concepts + '\n' + vcf_column_names_line + '\n'

    # Print to verbose VCF, high quality VCF, Monitoring VCF and ID VCF files
    OUT_vcf.write(final_header)
    OUT_vcf_hq.write(final_header)
    VCF_monitoring.write(final_header)
    VCF_id.write(final_header)


# Write column names to tsv outfile
def print_tsv_header():
    tsv_header = ['CHROM', 'POS', 'REF', 'ALT', 'DP_HQ', 'REF_COUNT', 'ALT_COUNT', 'AF', 'Multi_REF', 'Multi_ALT',
                  'Multi_AF', 'Pval', 'FDR', 'STRAND', 'FISHER', 'OAAC', 'MAL',
                  'Sequence_Context', 'Homopolymer', 'FILTER']
    tsv_header = '\t'.join(tsv_header) + '\n'
    OUT_tsv.write(tsv_header)
    OUT_tsv_hq.write(tsv_header)
    if bool(monitoring_variants):
        TSV_monitoring.write(tsv_header)
        TSV_id.write(tsv_header)


# Main method ---------------------------------------------------------------------------------------------

# Parse arguments
parser = argparse.ArgumentParser(description='Getting barcodes in fastq file and labels')
parser.add_argument('-i', '--infile', type=str, help='Table in tsv format (stats.txt)', required=True)
parser.add_argument('-s', '--single', type=str, help='Table in tsv format (dedup_DP1.tsv)', required=True)
parser.add_argument('-ref', '--reference', type=str, help='Reference fasta file', required=True)
parser.add_argument('-o', '--outfile', type=str, help='Output file in VCF format', required=True)
parser.add_argument('-m', '--monitoring', type=str, default='', help='Variants for monitoring or IDing in BED format')
parser.add_argument('-tID', '--tumorid', type=str, default='Tumor', help='Tumor sample id')
parser.add_argument('-cov', '--min_COV', type=int, default=10, help='Minimum Coverage')
parser.add_argument('-ac', '--min_AC', type=int, default=3, help='Minimum reads supporting alternative allele')
parser.add_argument('-af', '--min_AF', type=float, default=0.001, help='Minimum alternative allele frequency')
parser.add_argument('-dist', '--min_DIST', type=int, default=20, help='Minimum distance between variants')
parser.add_argument('-sb', '--strand_bias', type=int, choices=[0, 1], default=1, help='Fisher strand bias test')
parser.add_argument('-mrd', '--mrd', type=int, choices=[0, 1], default=1, help='Print Minimal Residual Disease')
parser.add_argument('-tmpdir', '--tmpdir', default=None, help='Folder for temp files')

# Store arguments
args = parser.parse_args()
input_file = args.infile
singleton_umi_file = args.single
sample_id = args.tumorid
reference_file = args.reference
monitoring_tsv = args.monitoring
min_COV = args.min_COV
min_AC = args.min_AC
min_AF = args.min_AF
min_DIST = args.min_DIST
strand = args.strand_bias


# Minimal Residual Disease containers
MRD_P = []  # List of pvalues at monitornig positions
MRD_DP = []  # List of reference counts at monitornig positions
MRD_ALT = []  # List of alternative counts at monitornig positions


# Read monitoring and ID variants (if file is provided by parameter)
monitoring_variants = dict()
if monitoring_tsv != '' and os.path.exists(monitoring_tsv):
    with open(monitoring_tsv) as fm:
        for current_line in fm:

            # Skip header lines in VCF files
            if current_line.startswith('##') or current_line.startswith('CHROM'):
                continue

            # extract fields of entry
            current_line = current_line.rstrip('\n')
            split_line = current_line.split("\t")
            current_chr = split_line[0]
            current_pos = split_line[1]
            current_type = split_line[2]
            current_ref = split_line[3]
            current_alt = split_line[4]
            variant_descriptors = [current_type, current_ref, current_alt]

            # Store entry in nested dictionary
            if current_chr not in monitoring_variants:
                monitoring_variants[current_chr] = dict()

            monitoring_variants[current_chr][current_pos] = variant_descriptors

# Open verbose output files
vcf_outfile_name = args.outfile
tsv_outfile_name = str(Path(vcf_outfile_name).with_suffix('.tsv'))
OUT_vcf = open(vcf_outfile_name, 'w')
OUT_tsv = open(tsv_outfile_name, 'w')

# Open high quality output files
file_base_folder = str(Path(vcf_outfile_name).parent)
file_base_name = str(Path(vcf_outfile_name).stem)
vcf_outfile_hq = file_base_folder + "/" + file_base_name + "_hq.vcf"
tsv_outfile_hq = file_base_folder + "/" + file_base_name + "_hq.tsv"
OUT_vcf_hq = open(vcf_outfile_hq, 'w')
OUT_tsv_hq = open(tsv_outfile_hq, 'w')

# Open Monitoring and ID variant files
TSV_monitoring = ''
TSV_id = ''
if bool(monitoring_variants):
    monitoring_tsv = file_base_folder + "/" + file_base_name + "_monitoring.tsv"
    id_tsv = file_base_folder + "/" + file_base_name + "_ID.tsv"
    TSV_monitoring = open(monitoring_tsv, 'w')
    TSV_id = open(id_tsv, 'w')

    monitoring_vcf = file_base_folder + "/" + file_base_name + "_monitoring.vcf"
    id_vcf = file_base_folder + "/" + file_base_name + "_ID.vcf"
    VCF_monitoring = open(monitoring_vcf, 'w')
    VCF_id = open(id_vcf, 'w')

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
OUT_vcf_hq.close()
OUT_tsv_hq.close()
if bool(monitoring_variants):
    TSV_monitoring.close()
    TSV_id.close()

exit(0)
