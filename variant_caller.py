# coding=utf-8

# Import Libraries
import argparse
import time
import numpy
import scipy.stats
from pathlib import Path
import pysam

numpy.seterr(divide='ignore')

def longestRun(s):
	if len(s) == 0: return 0
	runs = ''.join('*' if x == y else ' ' for x,y in zip(s,s[1:]))
	starStrings = runs.split()
	if len(starStrings) == 0: return 1
	return 1 + max(len(stars) for stars in starStrings)


def FrequentBase(s):
	L = list(s)
	
	# Get the most common element and the percentage respect the length of the sequence
	MAX = max([L.count(base) for base in set(L)])
	PERC = round(float(MAX)/len(s), 2)
	
	return(PERC)

def homopolymer_function(A):
	# A is the up or down-stream sequence
	if (A != '.'):
		# Get the longest k-mer
		Max_seq_l = longestRun(A)
		# Get the frequency of the more repeated base
		Max_base_freq = FrequentBase(A)
		
		if (Max_seq_l > 3 and Max_base_freq >= 0.8):
			FILTER = 1
		else:
			FILTER = 0
	else:
		FILTER = 0
		
	return(FILTER)


# Filter variants using various filter types
def variant_filter(variant, chrom, pos, mm, depth, depth_hq, dist, ref_t,up_context,down_context):

	# Declare variables
	gt, sig, alt_count, af, alt_count_p, alt_count_fdr, strand_counts, strand_bias, alt_count_other, out_adj = variant.split(":")

	# Reference and alternative base for this variant
	ref_base, alt_base = gt.split(">")

	# Get counts per strand
	alt_fwd, alt_rev, ref_fwd, ref_rev = strand_counts.split("-")

	# Homopolymer filter
	# Up-stream
	homopolymer_up = homopolymer_function(up_context)
	# Down-stream
	homopolymer_down = homopolymer_function(down_context)

	# Set filter based on various thresholds
	if float(alt_count) > 0 and mm > 0:
		filter_criteria = []

		# Filter: sequencing error FDR
		if float(alt_count_fdr) > 0.1:
			filter_criteria.append("Error")

		# Filter: homopolymer up-stream
		if homopolymer_up == 1:
			filter_criteria.append("LC_Upstream")

		# Filter: homopolymer down-stream
		if homopolymer_down == 1:
			filter_criteria.append("LC_Downstream")

		# Filter: minimum allele frequency
		if float(af) < min_AF:
			filter_criteria.append("Low_AF")

		# Filter: minimum alternative alleles per strand
		if (min(int(alt_fwd), int(alt_rev)) == 0 and
				min(int(ref_fwd), int(ref_rev)) != 0 and
				strand_counts == 1 and int(alt_count) > 3):
			filter_criteria.append("Strand_imbalanced")

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
			filter_criteria.append("Low_qual_pos")

		# Filter: contamination = high number of third or fourth allele
		if float(alt_count_other) / (float(mm)) > 0.3 or float(out_adj) < 0.1:
			filter_criteria.append("Variant_contamination")

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

	tsv = [alt_count, depth_hq, af, ref_t, alt_count_p, alt_count_fdr, strand_counts, strand_bias, alt_count_other, out_adj, homopolymer_up,homopolymer_down]

	variant_call = [ref_base, alt_base, tsv, concatenated_filter_string]

	return variant_call


# Parse aggregated pileup statistics (tsv) file with all regions of interest and call variants
def parse_pileup_statistics():

	# Open fasta file with pysam
	if reference_file != None:
		inFasta = pysam.FastaFile(reference_file)

	# Parse input file
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

			# Context sequences
			# Context regions (5 bases up- and down-stream)
			num_bases = 5
			if inFasta != None:
				try:
					context = inFasta.fetch(chrom, int(pos)-(num_bases+1), int(pos)+num_bases)
					context = context.upper()
					up_context = context[0:5]
					down_context = context[6:11]
				except:
					up_context = '.'
					down_context = '.'
			else:
				up_context = '.'
				down_context = '.'

			# Containers
			ref = []
			alt = []
			tsv = []

			# Apply variant filter method
			filter_string = []
			for variant in variant_candidates:
				variant_ref, variant_alt, variant_tsv, variant_filter_string = variant_filter(variant, chrom, pos, mm, dp, dp_hq, dist, ref_total,up_context,down_context)

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
						alt_count, dp_hq, ab_score, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, homopolymer_up, homopolymer_down = tsv[count]

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

						sample_index = [genotype, str(alt_count), str(dp_hq), str(ab_score), base_strand, str(fisher),
										str(alt_count_o), str(out_adj), str(alt_count_pad_j)]
						sample.append(':'.join(sample_index))

						p_val_list.append(float(alt_count_p))

					else:

						# Sample info
						alt_count, dp_hq, ab_score, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, homopolymer_up, homopolymer_down = tsv[count]

						alt_index = alt[count]
						alt_base.append(alt_index)
						alt_count_all.append(alt_count)

						if int(alt_count) > 0:
							genotype = alt_index
						else:
							genotype = reference_seq

						filter_field.append(filter_string[count])

						sample_index = [genotype, str(alt_count), str(dp_hq), str(ab_score),str(alt_count_pad_j),base_strand, str(fisher),
										str(alt_count_o), str(out_adj)]
						sample.append(':'.join(sample_index))

						p_val_list.append(float(alt_count_p))

				# Check numbers?
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
					format_field = ["GT", "Alt_Count", "DP", "AF", "Perror","Strand", "FS", "VCB", "Pvcb"]

					# Info column
					fields = ['Variant_dist=' + str(dist), 'Upstream='+str(up_context), 'Downstream='+str(down_context), 'PValue=' + str(alt_count_p)]

					# VCF variant line
					vcf_line = ['\t'.join(vcf_fields), filter_field[allele_idx], ';'.join(fields), ':'.join(format_field), sample[allele_idx]]
					vcf_line = '\t'.join(vcf_line)

					# TSV line
					tsv_line = [str(chrom), str(pos), str(reference_seq), str(alt_base[allele_idx]), str(dp_hq), str(ref_total),
								str(alt_count_all[allele_idx]), str(ab_score), alt_count_p, alt_count_pad_j, base_strand, fisher,
								alt_count_o, out_adj, str(up_context),str(homopolymer_up),str(down_context),str(homopolymer_down), filter_field[allele_idx]]
					tsv_line = '\t'.join(tsv_line)

					if int(alt_count) > 0:
						OUT_vcf.write(vcf_line + '\n')
						OUT_tsv.write(tsv_line + '\n')

			# Length of ref is <= 1
			else:

				# Common vcf columns
				reference_seq = ref[0]
				alt_base = alt[0]
				vcf_fields = [chrom, pos, snv_id, reference_seq, alt_base, qual]

				# Filter column
				filter_field = filter_string[0]

				# Sample info
				alt_count, dp_hq, ab_score, ref_total, alt_count_p, alt_count_pad_j, base_strand, fisher, alt_count_o, out_adj, homopolymer_up, homopolymer_down = tsv[0]

				# Info column
				fields = ['Variant_dist=' + str(dist), 'Upstream=' + str(up_context), 'Downstream=' + str(down_context),'PValue=' + str(alt_count_p)]

				# Format column
				format_field = ["GT", "Alt_Count", "DP", "AF", "Perror","Strand", "FS", "VCB", "Pvcb", ]

				# Append p-val for MRD calculation
				MRD.append(float(alt_count_p))

				# Get genotype
				if int(alt_count) > 0:
					genotype = alt_base
				else:
					genotype = reference_seq

				sample = [genotype, str(alt_count), str(dp_hq), str(ab_score),str(alt_count_pad_j), base_strand, str(fisher),
						  str(alt_count_o), str(out_adj)]

				# Compile VCF entry
				vcf_line = ['\t'.join(vcf_fields), filter_field, ';'.join(fields), ':'.join(format_field), ':'.join(sample)]
				vcf_line = '\t'.join(vcf_line)

				# Compile TSV entry
				tsv_line = [str(chrom), str(pos), str(ref[0]), str(alt[0]), str(dp_hq), str(ref_total),
							str(alt_count), str(ab_score), alt_count_p, alt_count_pad_j, base_strand, fisher,
							alt_count_o, out_adj, str(up_context),str(homopolymer_up),str(down_context),str(homopolymer_down), filter_field]
				tsv_line = '\t'.join(tsv_line)

				# Print variant candidate if at least one alternative read has been observed
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


# Print VCF header and column names
def print_vcf_header():
	# Define VCF header
	date = time.strftime("%d/%m/%Y")  # dd/mm/yyyy format
	vcf_format = "##fileformat=VCFv4.1"
	date_line = "##fileDate=%s" % date
	source = "##source=UKT-IMGAG UmiVar"
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
##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele read counts"
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternative allele fraction">
##FORMAT=<ID=Perror,Number=1,Type=Float,Description="Q-value to belong to the error rate distribution (beta binomial distribution) - After FDR correction">
##FORMAT=<ID=Strand,Number=2,Type=String,Description="Alleles in strands: Alternative forward, Alternative reverse, Reference forward, Reference reverse">
##FORMAT=<ID=FS,Number=1,Type=Float,Description="Fisher strand test (q-value)">
##FORMAT=<ID=VCB,Number=1,Type=Integer,Description="Variant Count bias, number of other different alternative alleles found">
##FORMAT=<ID=Pvcb,Number=1,Type=Float,Description="Q-value for other alternative alleles to belong to the error rate distribution (beta binomial distribution) - After FDR correction">"""

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
	tsv_header = ['CHROM', 'POS', 'REF', 'ALT', 'DP_HQ', 'REFt', 'ALT_COUNT', 'AF', 'P_VAL',
				  'P_VAL_adj', 'STRAND', 'FISHER', 'ALT_COUNT_o', 'P_VALo_adj', 'Up_context','Homopolymer_up','Down_context','Homopolymer_down', 'FILTER']
	tsv_header = '\t'.join(tsv_header)
	OUT_tsv.write(tsv_header + '\n')

# Main method ---------------------------------------------------------------------------------------------

# Parse arguments
parser = argparse.ArgumentParser(description='Getting barcodes in fastq file and labels')
parser.add_argument('-i', '--infile', type=str, help='Tsv table', required=True)
parser.add_argument('-tID', '--tumorid', type=str, default='Tumor', help='Tumor sample id', required=False)
parser.add_argument('-ref', '--reference', type=str, help='Reference fasta file which table was build', required=True)
parser.add_argument('-o', '--outfile', type=str, help='Vcf output file', required=True)
parser.add_argument('-cov', '--min_COV', type=int, default=10, help='Minimum Coverage', required=False)
parser.add_argument('-ac', '--min_AC', type=int, default=3, help='Minimum reads supporting alternative allele',
					required=False)
parser.add_argument('-variant_dist', '--min_DIST', type=int, default=20,
					help='Minimum distance allowed between variants (to avoid clustered errors)', required=False)
parser.add_argument('-str', '--strand', type=int, choices=[0, 1], default=1,
					help='Strand bias test (Fisher test). 0 for turn it off', required=False)
parser.add_argument('-af', '--min_AF', type=float, default=0, help='Minimum allele frequency allowed', required=False)
parser.add_argument('-mrd', '--mrd', type=int, choices=[0, 1], default=1,
					help='Print Minimal Residual Disease [Default = 1]', required=False)
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


# Open output files
vcf_outfile_name = args.outfile
tsv_outfile_name = str(Path(vcf_outfile_name).with_suffix('.tsv'))
OUT_vcf = open(vcf_outfile_name, 'w')
OUT_tsv = open(tsv_outfile_name, 'w')

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
