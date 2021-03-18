# coding=utf-8

# Import Libraries
import collections
import glob
import os
import sys
import argparse
import tempfile
import time
import pysam
import subprocess


# Take time for speed check
start_time = time.time()


# Create argument parser
argument_parser = argparse.ArgumentParser("umiVar - variant calling with unique molecular barcodes")

# Add arguments
argument_parser.add_argument('-tbam', '--tbam', required=True, help='Tumor bam file')
argument_parser.add_argument('-nbam', '--nbam', default='', required=False, help='Normal bam file')
argument_parser.add_argument('-b', '--bed', default='', help='Bed file of the targeted regions. O-based')
argument_parser.add_argument('-r', '--ref', required=True, help='Reference genome - fasta')
argument_parser.add_argument('-hom', '--homopolymer', type=str, help='File with homopolymer positions in reference genome', required=True)
argument_parser.add_argument('-o', '--out_file', default='./Sample.vcf', help='Out vcf file')
argument_parser.add_argument('-p', '--param', default='', help='Beta-binomial parameters table')
argument_parser.add_argument('-mq', '--mq', default=30, help='Minimum mapping quality')
argument_parser.add_argument('-bq', '--bq', default=20, help='Minimum base quality')
argument_parser.add_argument('-d', '--dist', default=5, help='Minimum distance allowed between variants')
argument_parser.add_argument('-ac', '--ac', default=3, help='Minimum number of reads supporting a variant')
argument_parser.add_argument('-ns', '--num_sites', default=-1, help='Number of sites to be analysed')
argument_parser.add_argument('-str', '--strand', default=0, choices=['0', '1'],
                             help='Strand filter activation. 0 for deactivating the filter. Default [0]')
argument_parser.add_argument('-t', '--temp_dir', default='.', help='Temporary directory')

# Retrieve arguments
arguments = argument_parser.parse_args()

# Use tumor BAM filename as ID
tumor_id = arguments.tbam.replace('.bam', '')

# Get directory of running script
script_directory = os.path.dirname(os.path.realpath(sys.argv[0]))

# Get output directory based on out_file (default is ./)
out_directory = os.path.dirname(arguments.out_file)

# If the output directory does not exist, create it with all subdirectories
if not os.path.exists(out_directory):
    os.makedirs(out_directory)

# Create a temporary working directory
temp_directory_path = arguments.temp_dir
temp_directory = tempfile.mkdtemp(dir=arguments.temp_dir, prefix='umiVar_tmp.')


# Split the BAM file into files with 1 or 2 or 3 or 4 and more barcode duplicates
def split_bam(infile, outprefix, dp_count_outfile):
    """Split BAM files..."""

    # output BAM files
    out1_fp = outprefix + '_DP1.bam'
    out2_fp = outprefix + '_DP2.bam'
    out3_fp = outprefix + '_DP3.bam'
    out4_fp = outprefix + '_DP4.bam'

    dp_counter = collections.Counter()

    # input BAM file
    with pysam.AlignmentFile(infile, 'rb') as infile, \
            pysam.AlignmentFile(out1_fp, 'wb', template=infile) as out1, \
            pysam.AlignmentFile(out2_fp, 'wb', template=infile) as out2, \
            pysam.AlignmentFile(out3_fp, 'wb', template=infile) as out3, \
            pysam.AlignmentFile(out4_fp, 'wb', template=infile) as out4:

        # iterate over BAM file
        for aln in infile:
            # read DP tag
            try:
                dp = aln.get_tag('DP')
            except KeyError:
                continue

            # add to counter
            dp_counter.update([dp])

            # split by DP
            if dp >= 4:
                out4.write(aln)
            elif dp == 3:
                out3.write(aln)
            elif dp == 2:
                out2.write(aln)
            elif dp <= 1:
                out1.write(aln)

    if dp_count_outfile:
        with open(dp_count_outfile, 'w') as f:
            f.write('dp\tcount\n')
            for (dp, count) in sorted(dp_counter.items()):
                f.write('{}\t{}\n'.format(dp, count))


# Split BAM file into 4 BAM files with 1 or 2 or 3 or 4 and more barcode duplicates
out_prefix = temp_directory + '/dedup'
frequency = temp_directory + '/dp_freq.tsv'
split_bam(infile=arguments.tbam, outprefix=out_prefix, dp_count_outfile=frequency)


# Debug speed
end_time = time.time()
time_taken = end_time - start_time
start_time = end_time
print('\nElapsed time for splitting BAM file: ' + str(time_taken) + "\n\n")


# Go over all of the split BAM files
for file in sorted(glob.glob(temp_directory + '/dedup*.bam')):

    # Index BAM files
    try:
        subprocess.run("samtools index " + file, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    # Create output file names for mpileup and pileup2tsv
    f_pileup = file.replace('.bam', '.pileup')
    f_tsv = file.replace('.bam', '.tsv')

    # Generate pileup and transform
    mpileup_command = 'samtools mpileup' + \
                      ' -d 0 ' + \
                      ' -f ' + arguments.ref + \
                      ' -Q 1 ' + \
                      ' -x ' + file + \
                      ' -o ' + f_pileup

    if arguments.bed != '':
        mpileup_command = 'samtools mpileup' + \
                          ' -d 0 ' + \
                          ' -f ' + arguments.ref + \
                          ' -l ' + arguments.bed + \
                          ' -Q 1 ' + \
                          ' -x ' + file + \
                          ' -o ' + f_pileup

    try:
        subprocess.run(mpileup_command, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    # Debug speed
    end_time = time.time()
    time_taken = end_time - start_time
    start_time = end_time
    print('\nElapsed time mpileup on ' + file + ' : ' + str(time_taken) + "\n")

    # Remove the temporary bam and bam.bai files
    os.remove(file)
    os.remove(file + '.bai')

    # Parse pileup file and create statistics for each genomic position in tsv format
    pileup_parser_command = ' python3 ' + script_directory + '/pileup_parser.py' + \
                            ' --pileup ' + f_pileup + \
                            ' --outfile ' + f_tsv + \
                            ' --minBQ ' + str(arguments.bq) + \
                            ' --minDP 10'

    try:
        subprocess.run(pileup_parser_command, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

    # Debug speed
    end_time = time.time()
    time_taken = end_time - start_time
    start_time = end_time
    print('Elapsed time pileup_parser.py on ' + f_pileup + " : " + str(time_taken) + "\n\n")


# Compute parameters of the negative binomial error distribution for each type of nucleotide-changes
f_params = None

# If no parameter files is given, estimate the optimal parameters of the beta-binomial distribution
if arguments.param == '':
    f_params = temp_directory + '/parameters.txt'

    beta_binomial_parameters_command = 'Rscript ' + script_directory + '/R_scripts/beta_binomial_parameters.R' + \
                                       ' -t1 ' + temp_directory + '/dedup_DP1.tsv' + \
                                       ' -t2 ' + temp_directory + '/dedup_DP2.tsv' + \
                                       ' -t3 ' + temp_directory + '/dedup_DP3.tsv' + \
                                       ' -t4 ' + temp_directory + '/dedup_DP4.tsv' + \
                                       ' -o ' + f_params

    try:
        subprocess.run(beta_binomial_parameters_command, shell=True)
    except subprocess.CalledProcessError as error:
        print(error)

# Else use the given parameters for the beta-binomial distribution
else:
    f_params = arguments.param


# Debug speed
end_time = time.time()
time_taken = end_time - start_time
start_time = end_time
print('\n\nElapsed time calculating beta parameters: ' + str(time_taken))


# Compute quality statistics for each genome position for variant calling
statistical_functions_command = 'Rscript ' + script_directory + '/R_scripts/statistical_functions.R' + \
                                ' -t1 ' + temp_directory + '/dedup_DP1.tsv' + \
                                ' -t2 ' + temp_directory + '/dedup_DP2.tsv' + \
                                ' -t3 ' + temp_directory + '/dedup_DP3.tsv' + \
                                ' -t4 ' + temp_directory + '/dedup_DP4.tsv' + \
                                ' -params ' + f_params + \
                                ' -num ' + str(arguments.num_sites) + \
                                ' -o ' + temp_directory + '/stats.tsv'

try:
    subprocess.run(statistical_functions_command, shell=True)
except subprocess.CalledProcessError as error:
    print(error)


# Debug speed
end_time = time.time()
time_taken = end_time - start_time
start_time = end_time
print('\n\nElapsed time for R statistics: ' + str(time_taken))


# Generate a VCF file
variant_caller_command = 'python3 ' + script_directory + '/variant_caller.py' + \
                         ' -i ' + temp_directory + '/stats.tsv' + \
                         ' -tID ' + tumor_id + \
                         ' -ref ' + arguments.ref + \
                         ' -hom ' + arguments.homopolymer + \
                         ' -o ' + arguments.out_file + \
                         ' -cov 10' + \
                         ' -ac ' + str(arguments.ac) + \
                         ' --strand ' + str(arguments.strand) + \
                         ' -variant_dist ' + str(arguments.dist) + \
                         ' -tmpdir ' + temp_directory

try:
    subprocess.run(variant_caller_command, shell=True)
except subprocess.CalledProcessError as error:
    print(error)


# Debug speed
end_time = time.time()
time_taken = end_time - start_time
start_time = end_time
print('\n\nElapsed time for variant_caller.py: ' + str(time_taken) + "\n\n")


# Plot error rates with R
plot_error_rate_command = 'Rscript ' + script_directory + '/R_scripts/plot_error_rates.R' + \
                          ' -bc1 ' + temp_directory + '/dedup_DP1.tsv' + \
                          ' -bc2 ' + temp_directory + '/dedup_DP2.tsv' + \
                          ' -bc3 ' + temp_directory + '/dedup_DP3.tsv' + \
                          ' -bc4 ' + temp_directory + '/dedup_DP4.tsv' + \
                          ' -out ' + temp_directory

try:
    subprocess.run(plot_error_rate_command, shell=True)
except subprocess.CalledProcessError as error:
    print(error)


# Debug speed
end_time = time.time()
time_taken = end_time - start_time
start_time = end_time
print('\n\nElapsed time for R plot error rates: ' + str(time_taken) + "\n\n")
