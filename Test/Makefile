all: clean CFDNA_01 CFDNA_02 select_monitoring post_filtering mrd_calculation_01

##Test case for CI
ci: clean download_test_reference CFDNA_02 select_monitoring post_filtering mrd_calculation_01

clean:
	rm -rf umiVar_CFDNA_01 umiVar_CFDNA_02 select_monitoring_01 select_monitoring_02 select_monitoring_03 postFilter01

download_reference:
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	gunzip -c hs37d5.fa.gz |  sed -r 's/>/>chr/g' > GRCh37.fa
	rm hs37d5.fa.gz
	samtools faidx GRCh37.fa

download_test_reference:
	wget http://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.7.fa.gz
	gunzip -kc Homo_sapiens.GRCh37.dna.chromosome.7.fa.gz | sed -r 's/>/>chr/g' > chr7.fa
	samtools faidx chr7.fa
	rm Homo_sapiens.GRCh37.dna.chromosome.7.fa.gz
	# chr for monitoring SNP selection
	wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
	gunzip -kc Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz | sed -r 's/>/>chr/g' > chr10.fa
	samtools faidx chr10.fa
	rm Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz

CFDNA_01:  run_CFDNA_01 check_CFDNA_01

run_CFDNA_01:
	python3 ../umiVar.py -tbam CFDNA_01.bam -b CFDNA_01.bed -r GRCh37.fa -o umiVar_CFDNA_01/ -m CFDNA_01_monitoring.vcf -t umiVar_CFDNA_01/tmp/

check_CFDNA_01:
	-diff -u CFDNA_01.tsv.out umiVar_CFDNA_01/CFDNA_01.tsv
	-diff -u CFDNA_01_hq.vcf.out umiVar_CFDNA_01/CFDNA_01_hq.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u CFDNA_01_monitoring.tsv.out umiVar_CFDNA_01/CFDNA_01_monitoring.tsv 
	-diff -u CFDNA_01_monitoring.vcf.out umiVar_CFDNA_01/CFDNA_01_monitoring.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u CFDNA_01_ID.tsv.out umiVar_CFDNA_01/CFDNA_01_ID.tsv 
	-diff -u CFDNA_01_ID.vcf.out umiVar_CFDNA_01/CFDNA_01_ID.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u CFDNA_01.mrd.out umiVar_CFDNA_01/CFDNA_01.mrd

CFDNA_02: run_CFDNA_02 check_CFDNA_02

run_CFDNA_02:
	python3 ../umiVar.py -tbam CFDNA_02.bam -b CFDNA_02.bed -r chr7.fa -o umiVar_CFDNA_02/ -m CFDNA_02_monitoring.vcf -t umiVar_CFDNA_02/tmp/

check_CFDNA_02:
	diff -u CFDNA_02.tsv.out umiVar_CFDNA_02/CFDNA_02.tsv
	diff -u CFDNA_02_hq.vcf.out umiVar_CFDNA_02/CFDNA_02_hq.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	diff -u CFDNA_02_monitoring.tsv.out umiVar_CFDNA_02/CFDNA_02_monitoring.tsv 
	diff -u CFDNA_02_monitoring.vcf.out umiVar_CFDNA_02/CFDNA_02_monitoring.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	diff -u CFDNA_02_ID.tsv.out umiVar_CFDNA_02/CFDNA_02_ID.tsv 
	diff -u CFDNA_02_ID.vcf.out umiVar_CFDNA_02/CFDNA_02_ID.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	diff -u CFDNA_02.mrd.out umiVar_CFDNA_02/CFDNA_02.mrd


select_monitoring: run_select_monitoring_01 run_select_monitoring_02 run_select_monitoring_03 check_select_monitoring_01 check_select_monitoring_02 check_select_monitoring_03

run_select_monitoring_01:
	python3 ../select_monitoring_variants.py -g selectMonitoring.GSvar -v selectMonitoring.vcf -r chr10.fa -o selectMonitoring01

check_select_monitoring_01:
	diff -u selectMonitoring_monitoring_01.bed.out selectMonitoring01/monitoring.bed
	diff -I "^##" -u selectMonitoring_monitoring_01.vcf.out selectMonitoring01/monitoring.vcf
	
run_select_monitoring_02:
	# ignore INDELs
	python3 ../select_monitoring_variants.py -g selectMonitoring.GSvar -v selectMonitoring.vcf -r chr10.fa -o selectMonitoring02 -i -n 50

check_select_monitoring_02:
	diff -u selectMonitoring_monitoring_02.tsv.out selectMonitoring02/monitoring.tsv
	diff -u selectMonitoring_ranked_02.tsv.out selectMonitoring02/ranked.tsv
	

run_select_monitoring_03:
	# ignore INDELs and off-target
	python3 ../select_monitoring_variants.py -g selectMonitoring.GSvar -v selectMonitoring.vcf -r chr10.fa -o selectMonitoring03 -i -t -n 50

check_select_monitoring_03:
	diff -u selectMonitoring_monitoring_03.tsv.out selectMonitoring03/monitoring.tsv
	diff -u selectMonitoring_ranked_03.tsv.out selectMonitoring03/ranked.tsv

post_filtering: run_post_filter_01 check_post_filter_01

run_post_filter_01:
	mkdir -p postFilter01
	python3 ../cfDNA_postfiltering.py CFDNA03_01.GSvar,CFDNA04_01.GSvar,CFDNA05_01.GSvar,CFDNA06_01.GSvar \
		postFilter01/CFDNA03_01.GSvar,postFilter01/CFDNA04_01.GSvar,postFilter01/CFDNA05_01.GSvar,postFilter01/CFDNA06_01.GSvar \
		--tumor_samples DNA01_01-DNA02_01.GSvar

check_post_filter_01:
	diff -u CFDNA03_01.GSvar.out postFilter01/CFDNA03_01.GSvar
	diff -u CFDNA03_01_combined.GSvar.out postFilter01/CFDNA03_01_combined.GSvar
	diff -u CFDNA04_01.GSvar.out postFilter01/CFDNA04_01.GSvar
	diff -u CFDNA05_01.GSvar.out postFilter01/CFDNA05_01.GSvar
	diff -u CFDNA06_01.GSvar.out postFilter01/CFDNA06_01.GSvar

mrd_calculation_01: run_mrd_calculation_01 check_mrd_calculation_01

run_mrd_calculation_01:
	mkdir -p test_output
	mkdir -p test_output/mrdCalculation01
	python3 ../calculateMRD.py data_input/mrdCalculation01 data_input/mrdCalculation01/cfdna_panel.vcf test_output/mrdCalculation01/mrdCalculation01_bg.mrd \
		test_output/mrdCalculation01/mrdCalculation01_monitoring_counts.tsv
	python3 ../calculateMRD.py data_input/mrdCalculation01 data_input/mrdCalculation01/cfdna_panel.vcf test_output/mrdCalculation01/mrdCalculation01_bg_unfiltered.mrd \
		test_output/mrdCalculation01/mrdCalculation01_monitoring_counts_unfiltered.tsv --max_af 1.0 --keep_gonosomes --keep_indels

check_mrd_calculation_01:
	diff -u data_output/mrdCalculation01/mrdCalculation01_bg.mrd test_output/mrdCalculation01/mrdCalculation01_bg.mrd
	diff -u data_output/mrdCalculation01/mrdCalculation01_monitoring_counts.tsv test_output/mrdCalculation01/mrdCalculation01_monitoring_counts.tsv
	diff -u data_output/mrdCalculation01/mrdCalculation01_bg_unfiltered.mrd test_output/mrdCalculation01/mrdCalculation01_bg_unfiltered.mrd
	diff -u data_output/mrdCalculation01/mrdCalculation01_monitoring_counts_unfiltered.tsv test_output/mrdCalculation01/mrdCalculation01_monitoring_counts_unfiltered.tsv



