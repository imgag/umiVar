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
    gt, sig, alt_count, af, alt_count_p, alt_count_fdr, strand_counts, strand_bias, alt_count_other, out_adj = variant.split(":")

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
        if float(alt_count_p) > 0.1:
            filter_criteria.append("pvalue")

        # Filter: FDR not significant
        if float(alt_count_fdr) > 0.1:
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
        if float(alt_count_other) / (float(mm)) > 0.3 or float(out_adj) < 0.1:
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

    tsv = [alt_count, depth_hq, af, ref_t, alt_count_p, alt_count_fdr, strand_counts, strand_bias, alt_count_other,
           out_adj, homopolymer]

    variant_call = [ref_base, alt_base, tsv, concatenated_filter_string]

    return variant_call


# Print variants of various quality levels to separate files
def print_variants(vcf_line, tsv_line):
    # Extract fields of tsv entry
    tsv_entries = tsv_line.split("\t")

    # Check if printing is required (at least one alternative read)
    tsv_chr = tsv_entries[0]
    tsv_pos = int(tsv_entries[1])
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

    # Check if ID or monitoring variant - always write to separate files even if zero alt_count
    if bool(monitoring_variants):
        if tsv_chr in monitoring_variants.keys() and tsv_pos in monitoring_variants[tsv_chr].keys():

            # Write either ID or monitoring variant
            if monitoring_variants[tsv_chr][tsv_pos] == 'ID':
                OUT_id.write(tsv_line + '\n')
            elif monitoring_variants[tsv_chr][tsv_pos] == 'M':
                OUT_monitoring.write(tsv_line + '\n')


# Parse aggregated pileup statistics (tsv) file with all regions of interest and call variants
def parse_pileup_statistics():
    # Open fasta file with pysam
    in_fasta = pysam.FastaFile(reference_file)

    # Parse tsv file with counts of singleton UMIs (dedup_DP1.tsv)
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

    # Parse input tsv file (stats.txt)
    with open(input_file) as f1:
        for line in f1:
            line = line.rstrip('\n')

            # Ignore header lines
            if line.startswith('##'):
                continue

            # Read column names of input file
            elif line.startswith('CHROM'):
                header_file = line
                form = header_file.split('\t')

                # Assign column names to indices
                chrom_index = [i for i, x in enumerate(form) if x == "CHROM"][0]
                pos_index = [i for i, x in enumerate(form) if x == "POS"][0]
                dp_index = [i for i, x in enumerate(form) if x == "DP"][0]
                dp_hq_index = [i for i, x in enumerate(form) if x == "DP_HQ"][0]
                ref_fwd_index = [i for i, x in enumerate(form) if x == "REFf"][0]
                ref_rev_index = [i for i, x in enumerate(form) if x == "REFr"][0]
                dist_index = [i for i, x in enumerate(form) if x == "DIST"][0]
                mm_index = [i for i, x in enumerate(form) if x == "MM"][0]
                call_index = [i for i, x in enumerate(form) if x == "CALL"][0]

                continue

            # Parse rows with actual data
            fields = line.split("\t")

            # Assign values to common VCF columns
            chrom = fields[chrom_index]
            pos = fields[pos_index]
            ref_nucleotide = fields[2]
            snv_id = '.'
            ref_fwd = int(fields[ref_fwd_index])
            ref_rev = int(fields[ref_rev_index])
            ref_total = ref_fwd + ref_rev
            mm = int(fields[mm_index])
            dp = int(fields[dp_index])
            dp_hq = int(fields[dp_hq_index])
            qual = '.'
            dist = fields[dist_index]

            # Up to 3 possible alternative alleles per position
            variant_candidates = fields[call_index].split("|")

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

            # Containers
            ref = []
            alt = []
            tsv = []

            # Apply variant filter method
            filter_string = []
            for variant in variant_candidates:
                variant_ref, variant_alt, variant_tsv, variant_filter_string = variant_filter(variant, mm, dp, dp_hq,
                                                                                              dist, ref_total,
                                                                                              sequence_context)

                # Append variants
                ref.append(variant_ref)
                alt.append(variant_alt)
                tsv.append(variant_tsv)
                filter_string.append(variant_filter_string)

            # Why is that important?
            if len(ref) > 1:

                # Getting ref. For both SNP, Deletion and Insertion, the longer reference will represent the original reference
                reference_seq = max(ref, key=len)

                filter_field = []
                alt_base = []
                alt_count_all = []
                sample = []
                p_val_list = []
                n = len(ref)

                for count in range(0, n):
                    ref_index = ref[count]

                    if len(ref_index) != len(reference_seq):

                        # Sample info
                        alt_count, dp_hq, allele_frequency, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, homopolymer = tsv[count]

                        # Getting normalized alternative allele
                        alt_index = alt[count]
                        alt_temp = list(reference_seq)
                        alt_temp[0] = alt_index
                        alt_index = ''.join(alt_temp)
                        alt_base.append(alt_index)
                        alt_count_all.append(alt_count)

                        # Determine genotype:
                        if float(allele_frequency) > 0.8:
                            genotype = "1/1"
                        elif float(allele_frequency) > 0.3:
                            genotype = "0/1"
                        else:
                            genotype = "./."

                        filter_field.append(filter_string[count])

                        sample_index = [genotype, str(alt_count), str(dp_hq), str(allele_frequency), base_strand, str(fisher),
                                        str(alt_count_o), str(out_adj), str(alt_count_pad_j)]
                        sample.append(':'.join(sample_index))

                        p_val_list.append(float(alt_count_p))

                    else:

                        # Sample info
                        alt_count, dp_hq, allele_frequency, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, homopolymer = tsv[count]

                        alt_index = alt[count]
                        alt_base.append(alt_index)
                        alt_count_all.append(alt_count)

                        # Determine genotype:
                        if float(allele_frequency) > 0.8:
                            genotype = "1/1"
                        elif float(allele_frequency) > 0.3:
                            genotype = "0/1"
                        else:
                            genotype = "./."

                        filter_field.append(filter_string[count])

                        sample_index = [genotype, str(alt_count), str(dp_hq), str(allele_frequency), str(alt_count_pad_j),
                                        base_strand, str(fisher),
                                        str(alt_count_o), str(out_adj)]
                        sample.append(':'.join(sample_index))

                        p_val_list.append(float(alt_count_p))

                # Check numbers?
                if len(alt_base) != len(filter_field) or len(alt_base) != len(sample) or len(alt_base) != len(
                        alt_count_all):
                    raise ValueError("Number of alt alleles differs")

                # Combine p-values of this site
                p_merge = scipy.stats.combine_pvalues(p_val_list)[1]
                MRD.append(float(p_merge))

                # Split multi-allelic variants in single lines
                for allele_idx in range(len(alt_base)):
                    # Common vcf columns
                    vcf_fields = [chrom, pos, snv_id, reference_seq, alt_base[allele_idx], qual]

                    # Format column
                    format_field = ["GT", "Alt_Count", "DP", "AF", "Perror", "Strand", "FS", "VCB", "Pvcb"]

                    # Info column
                    fields = ['Variant_dist=' + str(dist), 'Vicinity=' + str(sequence_context),
                              'PValue=' + str(alt_count_p)]

                    # VCF variant line
                    vcf_line = ['\t'.join(vcf_fields), filter_field[allele_idx], ';'.join(fields),
                                ':'.join(format_field), sample[allele_idx]]
                    vcf_line = '\t'.join(vcf_line)

                    # Extract DP of singleton UMI fragments for alternative base at locus
                    locus_id = str(chrom) + "_" + str(pos)
                    alt_key = alt_base[allele_idx]
                    if len(reference_seq) > 1:
                        alt_key = 'D'
                    if len(alt_key) > 1:
                        alt_key = 'I'

                    singleton_ref_depth = singleton_dp[locus_id][ref_nucleotide]
                    singleton_alt_depth = singleton_dp[locus_id][alt_key]
                    multi_umi_ref_depth = int(ref_total) - int(singleton_ref_depth)
                    multi_umi_alt_depth = int(alt_count_all[allele_idx]) - int(singleton_alt_depth)
                    multi_umi_af = multi_umi_alt_depth / (multi_umi_ref_depth + multi_umi_alt_depth)

                    # TSV line
                    tsv_line = [str(chrom), str(pos), str(reference_seq), str(alt_base[allele_idx]), str(dp_hq),
                                str(ref_total), str(alt_count_all[allele_idx]), str(allele_frequency),
                                str(multi_umi_ref_depth), str(multi_umi_alt_depth), str(round(multi_umi_af, 6)),
                                alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj,
                                str(sequence_context), str(homopolymer), filter_field[allele_idx]]

                    tsv_line = '\t'.join(tsv_line)

                    # Print variant candidates
                    print_variants(vcf_line, tsv_line)

            # Length of ref is <= 1
            else:

                # Common vcf columns
                reference_seq = ref[0]
                alt_base = alt[0]
                vcf_fields = [chrom, pos, snv_id, reference_seq, alt_base, qual]

                # Filter column
                filter_field = filter_string[0]

                # Sample info
                alt_count, dp_hq, allele_frequency, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, homopolymer = \
                    tsv[0]

                # Info column
                fields = ['Variant_dist=' + str(dist), 'Vicinity=' + str(sequence_context),
                          'PValue=' + str(alt_count_p)]

                # Format column
                format_field = ["GT", "Alt_Count", "DP", "AF", "Perror", "Strand", "FS", "VCB", "Pvcb", ]

                # Append p-val for MRD calculation
                MRD.append(float(alt_count_p))

                # Determine genotype:
                if float(allele_frequency) > 0.8:
                    genotype = "1/1"
                elif float(allele_frequency) > 0.3:
                    genotype = "0/1"
                else:
                    genotype = "./."

                sample = [genotype, str(alt_count), str(dp_hq), str(allele_frequency), str(alt_count_pad_j), base_strand,
                          str(fisher),
                          str(alt_count_o), str(out_adj)]

                # Compile VCF entry
                vcf_line = ['\t'.join(vcf_fields), filter_field, ';'.join(fields), ':'.join(format_field), ':'.join(sample)]
                vcf_line = '\t'.join(vcf_line)

                # Extract DP of singleton UMI fragments for alternative base at locus
                locus_id = str(chrom) + "_" + str(pos)
                alt_key = alt[0]
                if len(ref[0]) > 1:
                    alt_key = 'D'
                if len(alt_key) > 1:
                    alt_key = 'I'

                singleton_ref_depth = singleton_dp[locus_id][ref_nucleotide]
                singleton_alt_depth = singleton_dp[locus_id][alt_key]
                multi_umi_ref_depth = int(ref_total) - int(singleton_ref_depth)
                multi_umi_alt_depth = int(alt_count) - int(singleton_alt_depth)
                multi_umi_af = multi_umi_alt_depth / (multi_umi_ref_depth + multi_umi_alt_depth)

                # Compile TSV entry
                tsv_line = [str(chrom), str(pos), str(ref[0]), str(alt[0]), str(dp_hq), str(ref_total),
                            str(alt_count), str(allele_frequency), str(multi_umi_ref_depth), str(multi_umi_alt_depth),
                            str(round(multi_umi_af, 6)), alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj,
                            str(sequence_context), str(homopolymer), filter_field]
                tsv_line = '\t'.join(tsv_line)

                # Print variant candidates
                print_variants(vcf_line, tsv_line)


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
##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele read counts">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternative allele fraction">
##FORMAT=<ID=Perror,Number=1,Type=Float,Description="Q-value to belong to the error rate distribution (beta binomial distribution) - After FDR correction">
##FORMAT=<ID=Strand,Number=1,Type=String,Description="Alleles in strands: Alternative forward, Alternative reverse, Reference forward, Reference reverse">
##FORMAT=<ID=FS,Number=1,Type=Float,Description="Fisher strand test (q-value)">
##FORMAT=<ID=VCB,Number=1,Type=Integer,Description="Variant Count bias, number of other different alternative alleles found">
##FORMAT=<ID=Pvcb,Number=1,Type=Float,Description="Q-value for other alternative alleles to belong to the error rate distribution (beta binomial distribution) - After FDR correction">"""

    # Write VCF header and column names
    vcf_column_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_id]
    vcf_column_names_line = "#" + '\t'.join(vcf_column_names)

    # Print to verbose VCF file
    OUT_vcf.write(vcf_format + '\n')
    OUT_vcf.write(date_line + '\n')
    OUT_vcf.write(source + '\n')
    OUT_vcf.write(reference + '\n')
    OUT_vcf.write(concepts + '\n')
    OUT_vcf.write(vcf_column_names_line + '\n')

    # Print to high quality VCF file
    OUT_vcf_hq.write(vcf_format + '\n')
    OUT_vcf_hq.write(date_line + '\n')
    OUT_vcf_hq.write(source + '\n')
    OUT_vcf_hq.write(reference + '\n')
    OUT_vcf_hq.write(concepts + '\n')
    OUT_vcf_hq.write(vcf_column_names_line + '\n')


# Write column names to tsv outfile
def print_tsv_header():
    tsv_header = ['CHROM', 'POS', 'REF', 'ALT', 'DP_HQ', 'REFt', 'ALT_COUNT', 'AF', 'Multi_REF', 'Multi_ALT',
                  'Multi_AF', 'P_VAL', 'P_VAL_adj', 'STRAND', 'FISHER', 'ALT_COUNT_o', 'P_VALo_adj',
                  'Sequence_Context', 'Homopolymer', 'FILTER']
    tsv_header = '\t'.join(tsv_header) + '\n'
    OUT_tsv.write(tsv_header)
    OUT_tsv_hq.write(tsv_header)
    if bool(monitoring_variants):
        OUT_monitoring.write(tsv_header)
        OUT_id.write(tsv_header)


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
monitoring_file = args.monitoring
min_COV = args.min_COV
min_AC = args.min_AC
min_AF = args.min_AF
min_DIST = args.min_DIST
strand = args.strand_bias

# Containers
MRD = []
monitoring_variants = dict()

# Read monitoring and ID variants (if file is provided by parameter)
if monitoring_file != '' and os.path.exists(monitoring_file):
    with open(monitoring_file) as fm:
        for current_line in fm:

            # extract fields of entry
            current_line = current_line.rstrip('\n')
            split_line = current_line.split("\t")
            current_chr = split_line[0]
            current_pos = int(split_line[1]) + 1
            current_type = split_line[3]

            # Store entry in nested dictionary
            if current_chr not in monitoring_variants:
                monitoring_variants[current_chr] = dict()

            monitoring_variants[current_chr][current_pos] = current_type

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
OUT_monitoring = ''
OUT_id = ''
if bool(monitoring_variants):
    monitoring_file = file_base_folder + "/" + file_base_name + "_monitoring.tsv"
    id_file = file_base_folder + "/" + file_base_name + "_ID.tsv"
    OUT_monitoring = open(monitoring_file, 'w')
    OUT_id = open(id_file, 'w')

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
    OUT_monitoring.close()
    OUT_id.close()

exit(0)
