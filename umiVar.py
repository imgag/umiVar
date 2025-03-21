# coding=utf-8

# Import Libraries
import collections
import glob
import os
import sys
import argparse
import time
import pysam
import subprocess


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
        with open(dp_count_outfile, 'w') as dp_file:
            dp_file.write('dp\tcount\n')
            for (dp, count) in sorted(dp_counter.items()):
                dp_file.write('{}\t{}\n'.format(dp, count))


# Parse custom config file
def parse_config_file(config_file_path):
    if not os.path.exists(config_file_path):
        # if no settings file found use system binaries
        print("No settings file found! Using system versions of python, R and samtools.")
        return {"python": "python3", "R": "Rscript", "samtools": "samtools"}
    config = dict()
    print("Settings file found. Using the following binaries:")
    with open(config_file_path) as config_file:
        for line in config_file:
            # skip empty or comment lines
            if line.startswith("#") or line.startswith("[") or (line.strip() == ""):
                continue
            parsed_line = line.split('=')
            config[parsed_line[0].strip()] = parsed_line[1].strip()
            print("\t" + line.strip())

    # check mandatory settings
    mandatory = ["python", "R", "samtools"]
    for setting in mandatory:
        if setting not in config.keys():
            raise ValueError("Required setting '" + setting + "' not found in config file!")
    return config


# Parse environment variables for binary paths
def parse_env_config(config):
    for key in config.keys():
        value = os.getenv("umiVar_" + key + "_binary")
        if not value is None:
            print("Binary path for " + key + " found in environment variables: " + value)
            config[key] = value
    return config


# Main function
def main():
    # Take time for speed check
    start_time = time.time()

    # Create argument parser
    argument_parser = argparse.ArgumentParser("umiVar - variant calling with unique molecular barcodes")

    # Add arguments
    argument_parser.add_argument('-tbam', '--tbam', required=True, help='Tumor bam file')
    argument_parser.add_argument('-nbam', '--nbam', default='', required=False, help='Normal bam file')
    argument_parser.add_argument('-r', '--ref', required=True, help='Reference genome - fasta')
    argument_parser.add_argument('-b', '--bed', default='', help='Bed file of the targeted regions. O-based')
    argument_parser.add_argument('-m', '--monitoring', default='',
                                 help='VCF file with genomic positions for monitoring or sample-IDing.')
    argument_parser.add_argument('-o', '--out_folder', default='',
                                 help='Output folder. Will be created if not existing')
    argument_parser.add_argument('-p', '--param', default='', help='Beta-binomial parameters table')
    argument_parser.add_argument('-mq', '--mq', default=30, help='Minimum mapping quality')
    argument_parser.add_argument('-bq', '--bq', default=20, help='Minimum base quality')
    argument_parser.add_argument('-d', '--dist', default=5, help='Minimum distance between variants')
    argument_parser.add_argument('-ac', '--ac', default=5, help='Minimum number of reads supporting a variant')
    argument_parser.add_argument('-af', '--af', default=0.001, help='Minimum fraction of reads supporting a variant')
    argument_parser.add_argument('-ns', '--num_sites', default=-1, help='Number of sites to be analysed')
    argument_parser.add_argument('-sb', '--strand_bias', default=0, choices=['0', '1'],
                                 help='Fisher strand bias filter. Default [0]')
    argument_parser.add_argument('-t', '--temp_dir', default='.', help='Temporary directory')
    argument_parser.add_argument('-kt', '--keep_temp', action='store_true', help='Don\'t delete temporary directory')

    # Retrieve arguments
    arguments = argument_parser.parse_args()

    # Check if mandatory input files exist
    if not os.path.exists(arguments.tbam):
        print("Input bam file (-tbam) does not exist.")
        exit(1)
    if not os.path.exists(arguments.ref):
        print("Reference genome fasta file (-ref) does not exist.")
        exit(1)

    # Use tumor BAM filename as ID
    tumor_id = os.path.basename(os.path.splitext(arguments.tbam)[0])

    # Get directory of running script
    script_directory = os.path.dirname(os.path.realpath(sys.argv[0]))

    # Define and create output directory
    out_dir = arguments.out_folder

    # binary paths
    config = parse_config_file(os.path.join(script_directory, "settings.ini"))
    config = parse_env_config(config)

    # No output folder specified: create folder in current working directory
    if not out_dir:
        out_dir = "umiVar_out_" + tumor_id

    if os.path.exists(out_dir):
        print("Output directory already exists. Files might get overwritten.")

    os.makedirs(out_dir, exist_ok=True)

    # Split BAM file into 4 BAM files with 1 or 2 or 3 or 4 and more barcode duplicates
    out_prefix = out_dir + '/dedup'
    frequency = out_dir + '/dp_freq.tsv'
    split_bam(infile=arguments.tbam, outprefix=out_prefix, dp_count_outfile=frequency)

    # Debug speed
    end_time = time.time()
    time_taken = end_time - start_time
    start_time = end_time
    print('\nElapsed time for splitting BAM file: ' + str(time_taken) + "\n\n")

    # Go over all of the split BAM files
    for file in sorted(glob.glob(out_dir + '/dedup*.bam')):

        # Index BAM files
        try:
            subprocess.run(config["samtools"] + " index " + file, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            print(error)
            exit(1)

        # Create output file names for mpileup and pileup2tsv
        f_pileup = file.replace('.bam', '.pileup')
        f_tsv = file.replace('.bam', '.tsv')

        # Generate pileup and transform
        mpileup_command = config["samtools"] + ' mpileup' + \
                          ' -d 0 ' + \
                          ' -f ' + arguments.ref + \
                          ' -Q 1 ' + \
                          ' -x ' + file + \
                          ' -o ' + f_pileup

        if arguments.bed != '':
            mpileup_command = config["samtools"] + ' mpileup' + \
                              ' -d 0 ' + \
                              ' -f ' + arguments.ref + \
                              ' -l ' + arguments.bed + \
                              ' -Q 1 ' + \
                              ' -x ' + file + \
                              ' -o ' + f_pileup

        try:
            subprocess.run(mpileup_command, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            print(error)
            exit(1)

        # Debug speed
        end_time = time.time()
        time_taken = end_time - start_time
        start_time = end_time
        print('\nElapsed time mpileup on ' + file + ' : ' + str(time_taken) + "\n")

        # Remove the temporary bam and bam.bai files
        os.remove(file)
        os.remove(file + '.bai')

        # Parse pileup file and create statistics for each genomic position in tsv format
        pileup_parser_command = config["python"] + " " + script_directory + '/pileup_parser.py' + \
                                ' --pileup ' + f_pileup + \
                                ' --outfile ' + f_tsv + \
                                ' --minBQ ' + str(arguments.bq) + \
                                ' --minDP 10'

        try:
            subprocess.run(pileup_parser_command, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            print(error)
            exit(1)

        # Debug speed
        end_time = time.time()
        time_taken = end_time - start_time
        start_time = end_time
        print('Elapsed time pileup_parser.py on ' + f_pileup + " : " + str(time_taken) + "\n\n")

    # Compute parameters of the beta-binomial error distribution for each type of nucleotide-changes
    # If no parameter files is given, estimate the optimal parameters of the beta-binomial distribution
    if arguments.param == '':
        f_params = out_dir + '/beta_binom_parameters.txt'

        beta_binomial_parameters_command = config["R"] + " " + \
                                           script_directory + '/R_scripts/beta_binomial_parameters.R' + \
                                           ' -t1 ' + out_dir + '/dedup_DP1.tsv' + \
                                           ' -t2 ' + out_dir + '/dedup_DP2.tsv' + \
                                           ' -t3 ' + out_dir + '/dedup_DP3.tsv' + \
                                           ' -t4 ' + out_dir + '/dedup_DP4.tsv' + \
                                           ' -o ' + f_params

        try:
            subprocess.run(beta_binomial_parameters_command, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            print(error)
            exit(1)

    # else use the given parameters for the beta-binomial distribution
    else:
        f_params = arguments.param

    # Debug speed
    end_time = time.time()
    time_taken = end_time - start_time
    start_time = end_time
    print('\n\nElapsed time calculating beta parameters: ' + str(time_taken))

    # Compute quality statistics for each genome position for variant calling
    statistical_functions_command = config["R"] + " " + script_directory + '/R_scripts/statistical_functions.R' + \
                                    ' -t1 ' + out_dir + '/dedup_DP1.tsv' + \
                                    ' -t2 ' + out_dir + '/dedup_DP2.tsv' + \
                                    ' -t3 ' + out_dir + '/dedup_DP3.tsv' + \
                                    ' -t4 ' + out_dir + '/dedup_DP4.tsv' + \
                                    ' -params ' + f_params + \
                                    ' -num ' + str(arguments.num_sites) + \
                                    ' -o ' + out_dir + '/stats.tsv'

    try:
        subprocess.run(statistical_functions_command, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    # Debug speed
    end_time = time.time()
    time_taken = end_time - start_time
    start_time = end_time
    print('\n\nElapsed time for R statistics: ' + str(time_taken))

    # Generate a VCF file
    vcf_file = out_dir + "/" + tumor_id + ".vcf"

    # Variant calling in complete ROI
    if arguments.monitoring == '':
        variant_caller_command = config["python"] + " " + script_directory + '/variant_caller.py' + \
                                 ' -i ' + out_dir + '/stats.tsv' + \
                                 ' -s ' + out_dir + '/dedup_DP1.tsv' + \
                                 ' -tID ' + tumor_id + \
                                 ' -ref ' + arguments.ref + \
                                 ' -o ' + vcf_file + \
                                 ' -cov 10' + \
                                 ' -ac ' + str(arguments.ac) + \
                                 ' -af ' + str(arguments.af) + \
                                 ' -sb ' + str(arguments.strand_bias) + \
                                 ' -dist ' + str(arguments.dist) + \
                                 ' -tmpdir ' + out_dir

    # Variant calling in complete ROI, with focus on SNVs for monitoring or IDing
    else:
        variant_caller_command = config["python"] + " " + script_directory + '/variant_caller.py' + \
                                 ' -i ' + out_dir + '/stats.tsv' + \
                                 ' -s ' + out_dir + '/dedup_DP1.tsv' + \
                                 ' -tID ' + tumor_id + \
                                 ' -ref ' + arguments.ref + \
                                 ' -m ' + arguments.monitoring + \
                                 ' -o ' + vcf_file + \
                                 ' -cov 10' + \
                                 ' -ac ' + str(arguments.ac) + \
                                 ' -af ' + str(arguments.af) + \
                                 ' -sb ' + str(arguments.strand_bias) + \
                                 ' -dist ' + str(arguments.dist) + \
                                 ' -tmpdir ' + out_dir

    try:
        subprocess.run(variant_caller_command, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    # Debug speed
    end_time = time.time()
    time_taken = end_time - start_time
    start_time = end_time
    print('\n\nElapsed time for variant_caller.py: ' + str(time_taken) + "\n\n")

    # Plot error rates with R
    plot_error_rate_command = config["R"] + " " + script_directory + '/R_scripts/plot_error_rates.R' + \
                              ' -bc1 ' + out_dir + '/dedup_DP1.tsv' + \
                              ' -bc2 ' + out_dir + '/dedup_DP2.tsv' + \
                              ' -bc3 ' + out_dir + '/dedup_DP3.tsv' + \
                              ' -bc4 ' + out_dir + '/dedup_DP4.tsv' + \
                              ' -out ' + out_dir

    try:
        subprocess.run(plot_error_rate_command, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    # Debug speed
    end_time = time.time()
    time_taken = end_time - start_time
    print('\n\nElapsed time for R plot error rates: ' + str(time_taken) + "\n\n")


# Run program
main()
