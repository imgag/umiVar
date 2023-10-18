# Import Libraries
import argparse
import sys
import os
import pysam
import gzip

# Identifies and stores SNVs ------------------------------------------------------------------------------------------
class MonitoringVariant:

    # Constructor
    def __init__(self, input_vcf_file, input_gsv_file, out_folder, min_depth, min_alt, min_af, no_indels, num_var, no_off_target):

        # Input files and output folder
        self.input_vcf_file = input_vcf_file
        self.input_gsv_file = input_gsv_file
        self.out_folder = out_folder

        # Quality thresholds
        self.min_depth = min_depth
        self.min_alt = min_alt
        self.min_af = min_af
        self.no_indels = no_indels
        self.num_var = num_var
        self.no_off_target = no_off_target

        # Dictionaries storing information from input files
        self.gsv_gene = dict()
        self.gsv_score = dict()
        self.gsv_impact = dict()
        self.vcf = dict()

        # Dictionaries storing output
        self.good_mutation_counter = 0
        self.scoring_TSV = dict()
        self.scoring_VCF = dict()
        self.off_target_TSV = dict()
        self.off_target_VCF = dict()

    # Parse variant information from gsv file -------------------------------------------------------------------------
    def read_gsv(self):

        # Parse gsv file line by line
        with open(self.input_gsv_file) as gsv_file:
            for gsv_line in gsv_file:
                gsv_line = gsv_line.rstrip('\n')

                if gsv_line[0] == '#':
                    continue

                # Extract columns: chr, start, end, ref, obs, tumor_af, tumor_dp, normal_af, normal_dp, filter,
                # quality, gene, variant_type, coding_and_splicing, OMIM, ClinVar, HGMD, RepeatMasker,
                # dbSNP, 1000g, gnomAD, gnomAD_hom_hemi, gnomAD_sub, ESP_sub, phyloP, Sift, PolyPhen, fathmm-MKL,
                # CADD, REVEL, MaxEntScan, GeneSplicer, dbscSNV, COSMIC, NGSD_som_c, NGSD_som_p, NGSD_hom, NGSD_het
                # classification, classification_comment, validation, comment, gene_info, CGI_id, CGI_driver_statement
                # CGI_gene_role, CGI_transcript, CGI_gene, CGI_consequence, ncg_oncogene, ncg_tsg
                gsv_column = gsv_line.split("\t")

                # Score driverness
                driver = 0
                if 'known' in gsv_column[-7]:
                    driver = 2
                elif 'tier 1' in gsv_column[-7]:
                    driver = 1

                # Score role in oncogenesis
                role = 0
                if 'LoF' in gsv_column[-6] or \
                        'Act' in gsv_column[-6] or \
                        'ambiguous' in gsv_column[-6] or \
                        ('1' in gsv_column[-2] and 'na' not in gsv_column[-2]) or \
                        ('1' in gsv_column[-1] and 'na' not in gsv_column[-1]):
                    role = 1

                if 'not protein-affecting' in gsv_column[-7]:
                    role = 0

                #  Store score in dictionary
                locus = gsv_column[0] + "_" + gsv_column[1]
                onco_score = driver + role
                self.gsv_score[locus] = onco_score

                # Extract and store gene name and VEP Impact
                self.gsv_impact[locus] = 'MODIFIER'
                gene_column = gsv_column[11].split(",")
                if bool(gene_column):
                    self.gsv_gene[locus] = gene_column[0]
                    self.gsv_impact[locus] = gsv_column[13]

    # Parse variant information from vcf file -------------------------------------------------------------------------
    def evaluate_variants(self, reference_fasta):

        # Open output files
        fh_vcf = open(self.out_folder + '/monitoring.vcf', 'w')
        fh_tsv = open(self.out_folder + '/monitoring.tsv', 'w')
        fh_bed = open(self.out_folder + '/monitoring.bed', 'w')
        fh_rnk = open(self.out_folder + '/ranked.tsv', 'w')

        # Write tsv header files
        fh_tsv.write("CHROM\tPOS\tREF\tALT\tDepth\tREF_COUNT\tALT_COUNT\tAF\tFilter\tImpact\tGene\tContext\tScore\n")
        fh_rnk.write("CHROM\tPOS\tREF\tALT\tDepth\tREF_COUNT\tALT_COUNT\tAF\tFilter\tImpact\tGene\tContext\tScore\n")

        # Open reference genome file with pysam
        in_fasta = pysam.FastaFile(reference_fasta)

        # Parse vcf file line by line
        gzipped = False
        gzip_magic_number = "1f8b"
        fh_vcf_in = open(self.input_vcf_file)

        if fh_vcf_in.read(2).encode("utf-8").hex() == gzip_magic_number:
            fh_vcf_in = gzip.open(self.input_vcf_file)
        else:
            fh_vcf_in = open(self.input_vcf_file)

        caller = "unknown";

        with fh_vcf_in as vcf_file:
            for vcf_line in vcf_file:
                vcf_line = vcf_line.rstrip('\n')

                # VCF Header lines are simply printed to output vcf file
                if vcf_line[0] == '#':
                    fh_vcf.write(vcf_line + '\n')
                    
                    #check for caller:
                    
                    if vcf_line.startswith("##source="):
                        if "strelka" in vcf_line: 
                            caller = "strelka"
                        if "dragen" in vcf_line.lower(): 
                            caller = "dragen"
                    
                    continue
                   
                if caller == "unknown":
                    raise ValueError("Unknown caller for the VCF file, couldn't find caller line of supported caller (Strelka2, Dragen). Expected ##source=strelka or ##source=Dragen_somatic_calling")
                
                # Get VCF fields of a variant entry
                vcf_column = vcf_line.split("\t")

                # Get fields of the INFO and FORMAT columns
                info_column = vcf_column[7].split(";")
                format_column = vcf_column[8].split(":")

                # Get fields of the sample columns
                normal_sample = vcf_column[9].split(":")
                tumor_sample = vcf_column[10].split(":")

                # Store Info Abbreviations and content in dictionary
                info_pairs = dict()
                for info_string in info_column:
                    info_fields = info_string.split("=")
                    if len(info_fields) > 1:
                        info_pairs[info_fields[0]] = info_fields[1]

                # VEP annotation
                vep_annotation = info_pairs['CSQ'].split("|")

                # Store Sample fields in dictionary
                normal_table = dict()
                tumor_table = dict()

                for i in range(len(format_column)):
                    normal_table[format_column[i]] = normal_sample[i]
                    tumor_table[format_column[i]] = tumor_sample[i]

                # Check if indel
                is_indel = 0
                if len(vcf_column[3]) > 1 or len(vcf_column[4]) > 1:
                    is_indel = 1

                # Skip INDELS
                if is_indel == 1 and self.no_indels:
                    continue

                # Get depth of coverage for tumor and normal sample
                normal_dp = int(normal_table['DP'])
                tumor_dp = int(tumor_table['DP'])
                
                if caller == "strelka":
                # Compute REF_COUNT, ALT_COUNT and AF for SNVs
                    if is_indel == 0:
                        ref_format = vcf_column[3] + "U"
                        alt_format = vcf_column[4] + "U"
                        ref_t1, ref_t2 = tumor_table[ref_format].split(",")
                        alt_t1, alt_t2 = tumor_table[alt_format].split(",")
                        if (int(ref_t1) + int(alt_t1)) > 0:
                            af = int(alt_t1) / (int(ref_t1) + int(alt_t1))
                        else:
                            af = 0

                    # Compute REF_COUNT, ALT_COUNT and AF for indels
                    else:
                        ref_t1, ref_t2 = tumor_table['TAR'].split(",")
                        alt_t1, alt_t2 = tumor_table['TIR'].split(",")
                        tumor_dp = int(ref_t1) + int(alt_t1)
                        if tumor_dp > 0:
                            af = int(alt_t1) / tumor_dp
                        else:
                            af = 0
                elif caller == "dragen":
                    af = tumor_table['AF']
                    #depth of reference allele and alternate allele
                    ref_t1, alt_t1 = tumor_table['AD'].split(",")

                    af = float(af)
                    ref_t1 = int(ref_t1)
                    alt_t1 = int(alt_t1)
                    
                else:
                    raise ValueError("Unknown caller when trying to determine ref_count, alt_count and af.")

                # Generate locus ID from chr and pos (note: deletions in GSvar are not pos - 1
                locus_vcf = vcf_column[0] + "_" + vcf_column[1]
                locus_gsv = locus_vcf
                if len(vcf_column[3]) > 1:
                    gsv_pos = int(vcf_column[1]) + 1
                    locus_gsv = vcf_column[0] + "_" + str(gsv_pos)

                # Get gene name
                gene = ""
                if vep_annotation[3] != "":
                    gene = vep_annotation[3]
                if locus_gsv in self.gsv_gene:
                    gene = self.gsv_gene[locus_gsv]

                # Get Impact
                impact = "OTHER"
                if "LOW" in info_pairs['CSQ']:
                    impact = "LOW"
                if locus_gsv in self.gsv_impact and "LOW" in self.gsv_impact[locus_gsv]:
                    impact = "LOW"
                if "MODERATE" in info_pairs['CSQ']:
                    impact = "MODERATE"
                if locus_gsv in self.gsv_impact and "MODERATE" in self.gsv_impact[locus_gsv]:
                    impact = "MODERATE"
                if "HIGH" in info_pairs['CSQ']:
                    impact = "HIGH"
                if locus_gsv in self.gsv_impact and "HIGH" in self.gsv_impact[locus_gsv]:
                    impact = "HIGH"

                # Check for homopolymers and low-complexity regions in num_bases upstream and downstream
                num_bases = 5
                region_start = max(int(vcf_column[1]) - (num_bases + 1), 1)
                region_end = min(int(vcf_column[1]) + num_bases, in_fasta.get_reference_length(vcf_column[0]))

                try:
                    sequence_context = in_fasta.fetch(vcf_column[0], region_start, region_end)
                    sequence_context = sequence_context.upper()
                except FileNotFoundError:
                    sequence_context = 'NNNNNNNNNNN'

                hom_len = self.longest_homopolymer(sequence_context)
                distinct_bases, base_bias = self.frequent_base(sequence_context)

                print(sequence_context + "\t" + str(hom_len) + "\t" + str(distinct_bases) + "\t" + str(base_bias))

                # Score mutations
                score = 4 * af
                if af > 0.51:  # Very high AF indicates a germline variant
                    score -= 1

                if impact == "HIGH":  # High impact such as LoF, splice-defect or frameshift
                    score += 2
                elif impact == "MODERATE":  # Moderate impact such as missense variants
                    score += 1

                if vep_annotation[3] == "":  # Not in a gene
                    score -= 1

                if tumor_dp >= self.min_depth and int(alt_t1) >= self.min_alt:  # Good coverage
                    score += 1

                if impact == "OTHER":  # Malus for MODIFIER
                    score -= 1

                if locus_gsv in self.gsv_score:  # Driverness and role in oncogenesis
                    score += self.gsv_score[locus_gsv]

                # # Melanoma super-genes NRAS, KRAS, BRAF, NF1, TERT, CDKN2A, TP53 (replace by white-list file)
                # if gene == "NRAS" or gene == "KRAS" or gene == "BRAF" or gene == "NF1" \
                #         or gene == "TERT" or gene == "CDKN2A" or gene == "TP53":
                #     score += 1
                #
                # if gene == "TERT":
                #     if vcf_column[1] == 1295373 or vcf_column[1] == 1295250 or "12952" in vcf_column[1]:
                #         score += 2

                # Penalize low-complexity and homopolymer regions
                if hom_len >= 5:
                    score -= 1
                if distinct_bases < 3:
                    score -= 1
                if base_bias >= 0.8:
                    score -= 1
                
                # add score to VCF info column
                if vcf_column[7].strip() == "" or vcf_column[7].strip() == ".":
                    vcf_column[7] = "MonitoringScore=" + str(round(score, 3))
                else:
                    vcf_column[7] += ";MonitoringScore=" + str(round(score, 3))
                vcf_line = "\t".join(vcf_column)

                # Store info for high quality mutations (on-target and off-target separately)
                candidate = str(vcf_column[0]) + "\t" + str(vcf_column[1]) + "\t" + str(vcf_column[3]) + "\t" \
                            + str(vcf_column[4]) + "\t" + str(tumor_dp) + "\t" + str(ref_t1) + "\t" \
                            + str(alt_t1) + "\t" + str(round(af, 4)) + "\t" + str(vcf_column[6]) + "\t" \
                            + str(impact) + "\t" + str(gene) + "\t" + sequence_context + "\t" + str(round(score, 3))

                if "PASS" in vcf_column[6] or "." in vcf_column[6]:
                    self.good_mutation_counter += 1
                    self.scoring_TSV[candidate] = score
                    self.scoring_VCF[candidate] = vcf_line
                elif vcf_column[6] == "off-target" and normal_dp > 50:
                    self.off_target_TSV[candidate] = score
                    self.off_target_VCF[candidate] = vcf_line

        # Print results
        counter = 0
        for key in sorted(self.scoring_TSV, key=self.scoring_TSV.get, reverse=True):
            fh_rnk.write(key + "\n")

            # Print 30 SNPs for monitoring
            if counter < self.num_var:
                fh_tsv.write(key + "\n")
                fh_vcf.write(self.scoring_VCF[key] + "\n")
                columns = key.split("\t")
                start = int(columns[1]) - 1
                end = int(columns[1])
                fh_bed.write(columns[0] + "\t" + str(start) + "\t" + str(end) + "\n")
                counter += 1

        if not self.no_off_target:
            # add off-target reads to fill-up list
            for key in sorted(self.off_target_TSV, key=self.off_target_TSV.get, reverse=True):
                fh_rnk.write(key + "\n")

                # Fill-up 30 SNPs for monitoring
                if counter < self.num_var:
                    fh_tsv.write(key + "\n")
                    fh_vcf.write(self.off_target_VCF[key] + "\n")
                    columns = key.split("\t")
                    start = int(columns[1]) - 1
                    end = int(columns[1])
                    fh_bed.write(columns[0] + "\t" + str(start) + "\t" + str(end) + "\n")
                    counter += 1

        # Close output files
        fh_bed.close()
        fh_vcf.close()
        fh_tsv.close()
        fh_rnk.close()

    # Returns longest homopolymer
    @staticmethod
    def longest_homopolymer(sequence):
        if len(sequence) == 0:
            return 0

        runs = ''.join('*' if x == y else ' ' for x, y in zip(sequence, sequence[1:]))
        star_strings = runs.split()
        if len(star_strings) == 0:
            return 1

        return 1 + max(len(stars) for stars in star_strings)

    # Returns number of distinct bases and most common base frequency in a string
    @staticmethod
    def frequent_base(sequence):
        nucleotide_list = list(sequence.upper())

        # Most common base frequency
        most_common_base = max([nucleotide_list.count(base) for base in set(nucleotide_list)])
        most_common_base_percentage = round(float(most_common_base) / len(sequence), 2)

        # Number of distinct bases
        distinct_bases = 0
        for base in ['A', 'C', 'G', 'T']:
            if nucleotide_list.count(base) > 0:
                distinct_bases += 1

        return distinct_bases, most_common_base_percentage


# Main Method
def main():
    # Read parameters
    parser = argparse.ArgumentParser(description='Select variants for monitoring cancer treatment via liquid biopsy')
    parser.add_argument('-v', '--vcf', type=str, required=True, help='VCF file with all variants of a patient')
    parser.add_argument('-g', '--gsv', type=str, required=True, help='GSvar file with all variants of a patient')
    parser.add_argument('-r', '--ref', type=str, required=True, help='Reference genome file, e.g. GRCh38.fasta')
    parser.add_argument('-o', '--out', type=str, default='', help='Output directory. Must not exist.')
    parser.add_argument('-d', '--min_depth', type=int, default=50, help='Minimum depth at variant site.')
    parser.add_argument('-a', '--min_alt', type=int, default=5, help='Minimum alternative base count of variant.')
    parser.add_argument('-f', '--min_af', type=float, default=0.1, help='Minimum alternative allele frequency of variant.')
    parser.add_argument('-i', '--no_indels', action='store_true', help='Do not select INDELS as monitoring variants.')
    parser.add_argument('-n', '--num_var', type=int, default=30, help='Number of monitoring variants which will be selected.')
    parser.add_argument('-t', '--no_off_target', action='store_true', help='Do not fill-up list with off-target variants.')

    try:
        args = parser.parse_args()
    except IOError as io:
        print(io)
        sys.exit('Error reading parameters.')

    input_vcf_file = args.vcf
    input_gsv_file = args.gsv
    input_ref_file = args.ref
    out_dir = args.out
    min_depth = args.min_depth
    min_alt = args.min_alt
    min_af = args.min_af
    no_indels = args.no_indels
    num_var = args.num_var
    no_off_target = args.no_off_target

    # Create output folder
    # No output folder specified: create folder in current working directory
    if out_dir == '':
        out_dir = os.getcwd()
        out_dir += "/monitoring_xxx"

    # Full path to output folder specified
    elif out_dir[0] == "/":
        # Do not overwrite existing folders
        if os.path.exists(out_dir):
            print("Output directory already exists. Please specify new output folder.")
            exit(0)

    #  Relative path to output folder specified
    elif out_dir[0:2] == "./":
        out_dir = os.getcwd() + "/" + out_dir[2:]

    # Specified output directory is just a name
    else:
        out_dir = os.getcwd() + "/" + out_dir

    # Create output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        print("Warning: output directory already exists and will be overwritten!")
        # exit(0)

    # Instantiate MonitoringVariant object and run the variant evaluation
    evaluator = MonitoringVariant(input_vcf_file, input_gsv_file, out_dir, min_depth, min_alt, min_af, no_indels, num_var, no_off_target)
    evaluator.read_gsv()
    evaluator.evaluate_variants(input_ref_file)

    print('Finished')


# Run tool
main()
