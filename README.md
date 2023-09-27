# umiVar

Detection of ultra-low fraction variants using unique molecular barcodes and ultra-deep sequencing.\
PCR copies of the same original DNA molecule are recognized based on their identical barcode and \
mapping position. umiVar first computes corrected consensus reads with adjusted quality scores. \
Next, specific error models for the different nucleotide changes and for different levels of barcode\
correction are computed. The 'barcode-correction-level' depends on the number of PCR copies (duplicate\
reads with identical barcodes) that have been used for generating a consensus read, with quality levels\
of 1 (no copy, no correction), 2 (2 copies form the consensus), 3 copies and finally 4 or more copies.

Variant calling is based on the beta-binomial error models in combination with other well-known quality\
features (allele frequency, low complexity sequence, strand bias etc.). Variants can be called at fractions\
as low as 1 in 10000 reads (0.01%). In addition, umiVar provides various plots for analysing the error\
correction efficiency, and the PCR copy number distribution.


## Installation
Simply clone the repository

## Dependencies
Required python packages:
  - samtools
  - pysam
  - numpy
  - scipy
  - networkx

Required R packages:
  - VGAM
  - argparse
  - data.table
  - ggplot2
  - grid
  - gridExtra
  - scales
  - seqinr
  - survcomp

## Homopolymer detection

```
usage: homopolymer_finder.py [-h] -r REF [-l LENGTH] [-o OUTLIER]

Find all homopolymers of length n in reference sequence (fasta).

optional arguments:
  -h, --help            show this help message and exit
  -r REF, --ref REF     Reference fasta file.
  -l LENGTH, --length LENGTH
                        Minimum length of homopolymers to report. Default = 6
  -o OUTLIER, --outlier OUTLIER
                        Maximum number of outliers (other nucleotides) in homopolymers. Default = 0
```

Example: Detect all homopolymers of at least 6 consecutive equal nucleotides or 7mers with at least 6 equal nucleotides and a wildcard.
```
python homopolymer_finder.py -r GRCh37.fa -l 7 -o 1 > GRCh37.homopolymer.txt
```

## Barcode correction
```
usage: barcode_correction.py [-h] --infile INFILE --outfile OUTFILE [--barcodes {START,END,BOTH}] [--minBQ MINBQ] [--barcode_error BARCODE_ERROR] [--n]

Correcting BAM files using barcodes info

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       Input BAM file.
  --outfile OUTFILE     Output BAM file.
  --barcodes {START,END,BOTH}
                        Barcode position: START = 5' barcode; END = 3' barcode; BOTH = 5' and 3' barcodes. Default = BOTH
  --minBQ MINBQ         Minimum base quality to be considered. Default = 30
  --barcode_error BARCODE_ERROR
                        Maximum number of sequencing errors allowed in barcode sequence. Default = 0
  --n                   Use Ns instead of reducing base quality.
```
\
Example:
```
python barcode_correction.py --infile raw.bam --outfile barcode_corrected.bam --barcodes BOTH
```

## umiVar variant caller
```
usage: umiVar - variant calling with unique molecular barcodes [-h] -tbam TBAM [-nbam NBAM] [-b BED] -r REF -hom HOMOPOLYMER [-o OUT_FILE] [-p PARAM] [-mq MQ] [-bq BQ] [-d DIST]
                                                               [-ac AC] [-ns NUM_SITES] [-str {0,1}] [-t TEMP_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -tbam TBAM, --tbam TBAM
                        Tumor bam file
  -nbam NBAM, --nbam NBAM
                        Normal bam file
  -b BED, --bed BED     Bed file of the targeted regions. O-based
  -r REF, --ref REF     Reference genome - fasta
  -hom HOMOPOLYMER, --homopolymer HOMOPOLYMER
                        File with homopolymer positions in reference genome
  -o OUT_FILE, --out_file OUT_FILE
                        Out vcf file
  -p PARAM, --param PARAM
                        Beta-binomial parameters table
  -mq MQ, --mq MQ       Minimum mapping quality
  -bq BQ, --bq BQ       Minimum base quality
  -d DIST, --dist DIST  Minimum distance allowed between variants
  -ac AC, --ac AC       Minimum number of reads supporting a variant
  -ns NUM_SITES, --num_sites NUM_SITES
                        Number of sites to be analysed
  -str {0,1}, --strand {0,1}
                        Strand filter activation. 0 for deactivating the filter. Default [0]
  -t TEMP_DIR, --temp_dir TEMP_DIR
                        Temporary directory
```
\
Example:
```
python umiVar.py -tbam cfDNA.bam -b region.bed -r GRCh37.fa  -p beta_params.txt
```

### Settings file
You can also provide custom binaries for `python`, `R` and `samtools`. For that a settings file called `settings.ini`has to be created in the main directory of umiVar. If no `settings.ini`is available, umiVar uses the system-wide installed binaries.
\
Example content (also shown in `settings.default`):
```
python = [PATH_TO_PYTHON]/python3
R = [PATH_TO_R]/RScript
samtools = [PATH_TO_SAMTOOLS]/samtools
```
Additionally, these paths can also be provided by environment variables:
```
$umiVar_python_binary="[PATH_TO_PYTHON]/python3"
$umiVar_R_binary="[PATH_TO_R]/RScript"
$umiVar_samtools_binary="[PATH_TO_SAMTOOLS]/samtools"
```

