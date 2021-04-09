all: clean run check

clean:
	rm -rf cfdna_tmp.* test_out.* umiVar_tmp.*

run:
	python3 ../umiVar.py -tbam CFDNA_01.bam -b target.bed -r GRCh37.fa -o test_out_new.vcf -p p.txt -t .

check:
	-diff -u test_out_new.vcf test_evaluation/out.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u test_out_new.tsv test_evaluation/out.tsv
	-diff -u test_out_new.mrd test_evaluation/out.mrd
