all: clean run_personalized check_personalized

download_reference:
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	gunzip -c hs37d5.fa.gz |  sed -r 's/>/>chr/g' > GRCh37.fa
	rm hs37d5.fa.gz
	samtools faidx GRCh37.fa

clean:
	rm -rf umiVar_personalized

run_personalized:
	python3 ../umiVar.py -tbam CFDNA_01.bam -b CFDNA_01.bed -r GRCh37.fa -o umiVar_personalized/ -m CFDNA_01_monitoring.vcf -t umiVar_personalized/tmp/ -crb ../R-4.0.5/bin/Rscript

check_personalized:
	-diff -u CFDNA_01.tsv.out umiVar_personalized/CFDNA_01.tsv
	-diff -u CFDNA_01_hq.vcf.out umiVar_personalized/CFDNA_01_hq.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u CFDNA_01_monitoring.tsv.out umiVar_personalized/CFDNA_01_monitoring.tsv 
	-diff -u CFDNA_01_monitoring.vcf.out umiVar_personalized/CFDNA_01_monitoring.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u CFDNA_01_ID.tsv.out umiVar_personalized/CFDNA_01_ID.tsv 
	-diff -u CFDNA_01_ID.vcf.out umiVar_personalized/CFDNA_01_ID.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u CFDNA_01.mrd.out umiVar_personalized/CFDNA_01.mrd
