##mummer compare, check the mapping info of chr3 in M.polymorpha
cd /data/home/zhangfan/compare_genomoe_medicago/re_anno/Mpoly/
seqkit grep -f chr_list.txt Mpoly.fa > Mpoly.chr.fa ##only extract chromosome info for download analysis
python ../replace_chr_name.py chr_replace.txt Mpoly.chr.fa Mpoly.chr.use.fa ##match the chromosome number to MtA17.
seqkit grep -f chr_list.txt Mpoly.chr.use.fa > Mpoly.chr3.fa ##extract chr3 info, chr_list.txt contain:Chr3
nucmer --mum --t 64 -p Mpoly_MtA17 Mpoly.chr3.fa ../MtrunA17_5.0_chr.genome.fasta
delta-filter -1 -i 85 -l 5000 -q Mpoly_MtA17.delta > Mpoly_MtA17.delta-filter
show-coords -THrd Mpoly_MtA17.delta-filter >  Mpoly_MtA17.filtered.coords
##extract the gap info of chr3
python getgaps.py Mpoly.chr3.fa > gaps.gff3
##overlapping the Mpoly_MtA17.filtered.coords and gaps.gff3 info
##find the breakpoint
Chr3	.	gap	51900483	51900582	.	.	.	Name=gap9;size=100
##seperate the Chr3 into Chr3 and Chr7
echo -e "Chr3\t0\t51900483" | bedtools getfasta -fi Mpoly.chr3.fa -bed - -fo Mpoly.chr7.fa
sed -i 's/Chr3:0-51900483/Chr7/g' Mpoly.chr7.fa
echo -e "Chr3\t51900582\t93525394" | bedtools getfasta -fi Mpoly.chr3.fa -bed - -fo Mpoly.chr3_use.fa
sed -i 's/Chr3:51900582-93525394/Chr3/g' Mpoly.chr3_use.fa

seqkit grep -f chr_list.txt Mpoly.chr.use.fa > Mpoly.chr8.fa ##Edit chr_list.txt to reflect "Chr8"
seqkit seq -r -p --seq-type DNA Mpoly.chr8.fa > Mpoly.chr8_reverse.fa ##reverse chromosome, can be used for multiple chromosomes

##combine all chromosome
seqkit grep -f chr_list.txt Mpoly.chr.use.fa > Mpoly.chr_rest.fa ##extract the rest chromosomes
cat Mpoly.chr_rest.fa Mpoly.chr3_use.fa Mpoly.chr7.fa Mpoly.chr8_reversed.fa > Mpoly_merged_chromosomes.fa ##combine all chr
seqkit sort -n Mpoly_merged_chromosomes.fa -o Mpoly_merged_chromosomes_sort.fa ##sort chr

##we used the similar pipeline to reverse the chr for the rest genomes
