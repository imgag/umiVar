all: clean CFDNA_01 CFDNA_02

##Test case for CI
ci: clean download_test_reference CFDNA_02

clean:
	rm -rf umiVar_CFDNA_01 umiVar_CFDNA_02

download_reference:
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	gunzip -c hs37d5.fa.gz |  sed -r 's/>/>chr/g' > GRCh37.fa
	rm hs37d5.fa.gz
	samtools faidx GRCh37.fa

download_test_reference:
	wget http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.7.fa.gz
	gunzip -kc Homo_sapiens.GRCh37.dna.chromosome.7.fa.gz | sed -r 's/>/>chr/g' > chr7.fa
	samtools faidx chr7.fa
	rm Homo_sapiens.GRCh37.dna.chromosome.7.fa.gz

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

