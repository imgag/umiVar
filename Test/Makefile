all: clean run new check

clean:
	rm -rf cfdna_tmp.* test_out.*

run:
	python3 ../umiVar.py -tbam CFDNA_01.bam -b target.bed -r GRCh37.fa -o test_out_new.vcf -p p.txt -t .

check:
	-diff -u test_out_new.vcf out.vcf --ignore-matching-lines="^##fileDate" --ignore-matching-lines="^##reference"
	-diff -u test_out_new.tsv out.tsv
	-diff -u test_out_new.mrd out.mrd
