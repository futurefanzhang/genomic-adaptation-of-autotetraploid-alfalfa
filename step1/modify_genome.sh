##mummer compare, check the mapping info of chr3 in M.polymorpha
cd /data/home/zhangfan/compare_genomoe_medicago/re_anno/Mpoly/
nucmer --mum --t 64 -p Mpoly_MtA17 Mpoly.chr3.fa ../MtrunA17_5.0_chr.genome.fasta
delta-filter -1 -i 85 -l 5000 -q Mpoly_MtA17.delta > Mpoly_MtA17.delta-filter
show-coords -THrd Mpoly_MtA17.delta-filter >  Mpoly_MtA17.filtered.coords
##extract the gap info of chr3
python getgaps.py Mpoly.chr3.fa > gaps.gff3
##overlapping the Mpoly_MtA17.filtered.coords and gaps.gff3 info
##find the breakpoint
Chr3	.	gap	51900483	51900582	.	.	.	Name=gap9;size=100
