
Homework 3

Ria Deshpande  
November 9, 2020

### Summarize Genome Assembly
To summarize the genome assembly of the latest version of the Drosophila melanogaster genome, we first need to download the current file of the unaligned sequence as well as the md5sum file to verify integrity.

``` bash
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/md5sum.txt
md5sum -c md5sum.txt
```
Next, we need to identify the size of each of the sequences in the file using the tool faSize. Next, we can find the total number of nucleotides, find the total number of Ns, and the total number of sequences. These can be done as follows:
``` bash
all_url="ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz"
wget $all_url
faSize -detailed ~/HW3/dmel-all-chromosome-r6.36.fasta.gz > ~/HW3/all-chromosomes-namesizes.txt
sed '1,7!d' ~/HW3/all-chromosomes-namesizes.txt > all_chromosomes_sizes.txt #saves the first 7 lines in a different file
awk '{print $2}' all_chromosomes_sizes.txt > only_chromosomes_sizes.txt #prints only the second column, giving length of each chromosome in a new file
grep . ~/HW3/only_chromosomes_sizes.txt | paste -sd+ | bc #sums the column, which contains only length of each chromosome
grep -o 'N' ~/HW3/dmel-all-chromosome-r6.36.fasta.gz | wc -l #counts the number of N's in the whole file
wc -l all-chromosomes-namesizes.txt #counts the number of lines, assuming that each line represents one sequence
```


### Summarize Annotation Assembly
To summarize the annotation assembly, we use a series of small pipelines. Similar to above, we need to download the appropriate data and verify the integrity of the files. The answers are given in the comments below.

``` bash
wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/dmel-all-r6.36.gtf.gz
wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/md5sum.txt
md5sum -c md5sum.txt
gzip -d dmel-all-r6.36.gtf.gz
unzipGtfFile=dmel-all-r6.36.gtf

cat $unzipGtfFile | grep -v '^#' | cut -f3 | sort | uniq -c
# 33629 3UTR
# 46664 5UTR
# 162578 CDS
# 189268 exon
# 17875 gene
# 485 miRNA
# 30728 mRNA
# 3047 ncRNA
# 262 pre_miRNA
# 366 pseudogene
# 115 rRNA
# 300 snoRNA
# 32 snRNA
# 30812 start_codon
# 30754 stop_codon
# 312 tRNA

grep -v '^#' $unzipGtfFile |
awk -F '\t' '$1=="X" && $3=="gene"'  | wc -l
#2691

grep -v '^#' $unzipGtfFile |
awk -F '\t' '$1=="Y" && $3=="gene"'  | wc -l
#113

grep -v '^#' $unzipGtfFile |
awk -F '\t' '$1=="2L" && $3=="gene"'  | wc -l
#3516

grep -v '^#' $unzipGtfFile |
awk -F '\t' '$1=="2R" && $3=="gene"'  | wc -l
#3653

grep -v '^#' $unzipGtfFile |
awk -F '\t' '$1=="3L" && $3=="gene"'  | wc -l
#3486

grep -v '^#' $unzipGtfFile |
awk -F '\t' '$1=="3R" && $3=="gene"'  | wc -l
#4225

grep -v '^#' $unzipGtfFile |
awk -F '\t' '$1=="4" && $3=="gene"'  | wc -l
#114
```
