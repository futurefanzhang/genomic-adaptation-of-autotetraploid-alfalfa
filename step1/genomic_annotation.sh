###this part include genome annotation steps.
##we used WGAP (https://github.com/unavailable-2374/Genome-Wide-Annotation-Pipeline) of Shuo Cao's pipeline to do our annotation.
##download Medicago truncatula A17 5.0 annotation info (MtrunA17r5.0-ANR-EGN-r1.9.prot.fasta and Mt5-families.RepeatModeler-2.0.1.20200509.fa link: https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/)
micromamba activate WGAP
perl /data/home/zhangfan/compare_genomoe_medicago/Genome-Wide-Annotation-Pipeline-main/bin/GWAP.pl \
--genome Medicago_arabica.fa -1 ./ERR6688729/ERR6688729_1.fastq.gz -2 ERR6688729/ERR6688729_2.fastq.gz \
--protein /data/home/zhangfan/compare_genomoe_medicago/MtrunA17r5.0-ANR-EGN-r1.9.prot.fasta \
--out_prefix out_put --cpu 30 --gene_prefix Mara --Pfam_db /data/home/zhangfan/compare_genomoe_medicago/Genome-Wide-Annotation-Pipeline-main/Pfam-AB.hmm \
--augustus_species Mara --RM_lib /data/home/zhangfan/compare_genomoe_medicago/Medicago_arabica/Mt5-families.RepeatModeler-2.0.1.20200509.fa

##filter results
cat out_put.GeneModels.gff3 out_put.GeneModels_lowQuality.gff3 > combine_all_gene.gff3 ##把基因注释的结果合并(包含高质量和低质量)
GFF3Clear --genome /data/home/zhangfan/compare_genomoe_medicago/Medicago_lupulina/Medicago_lupulina.fa combine_all_gene.gff3 --gene_prefix Mlup > output_all_gene.gff3 #合并成新的gff3文件，基因的名称是Mlup开头
awk '$3 == "gene"' output_all_gene.gff3 > gene_only.gff3 #提取基因对应的行
awk '{ if (a[$4]) {a[$4]=a[$4]"\n"$0} else {a[$4]=$0} } END {for (i in a) if (a[i] ~ /\n/) print a[i]}' gene_only.gff3 > repeat_gene.gff3 #把起始位置重复的基因都挑选出来
perl /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/select_gene.pl ##运行这个脚步，删除的基因会保存在removed_genes.gff3中
perl extract_ID.pl ##用这个脚本把ID信息导出到remove_gene_list.txt中
grep -v -f remove_gene_list.txt output_all_gene.gff3 > filtered_dup_ID_all_gene.gff3 #去掉gff3文件中的重复基因，remove_gene_list.txt 为基因名称，每行一个ID对应的基因名
for i in $(cat species.txt);do gffread $(basename ${i}).gff3 -g $(basename ${i}).fa -x $(basename ${i}).cds.fa -y $(basename ${i}).pep.fa -w $(basename ${i}).trans.fa >> recode.txt 2>&1; done ##利用gff3和fa文件生成pep文件
###如果出现CDS 的end<start错误，把对应的ID的基因信息也从gff3文件去掉，然后重新做上面生成pep的一步
perl remove_error_gene.pl #利用这个脚步把recode.txt的基因导出到remove_gene_list2.txt
cat remove_gene_list.txt remove_gene_list2.txt > remove_gene_list3.txt
grep -v -f remove_gene_list3.txt output_all_gene.gff3 > filtered_dup_ID_all_gene.gff3 #去掉错误的基因名称
perl -ne 'if(/^>(.*)/){$gene=$1} else {if(/\./){print "$gene\n"}}' Medicago_lupulina.pep.fa > test.out.txt #利用perl命令查找包含.的基因名称
grep -v -f test.out.txt Medicago_lupulina.gff3 > final_gene.gff3


###Optional steps 
##去掉注释效果不好的基因，主要是Augustus_transcriptSupport_percentage=0;RNASeq_exon_base_depth_median=0和Augustus_intronSupport=0/0 这种类型，因为转录本支持比例是0，外显子深度是0，内含子也是0
awk '$3=="gene"{print $9}' ZM4_V2.gff3 | grep "Augustus_transcriptSupport_percentage=0" | grep "RNASeq_exon_base_depth_median=0"| grep "Augustus_intronSupport=0/0" |sed 's:.*;ID=Msa:Msa:' |  sed 's:\;.*::' | uniq > no_evidence.id.txt ##筛选出ID信息
grep -v -f no_evidence.id.txt ZM4_V2.gff3 > ZM4_V2_2.gff3 #筛选掉支持度低的基因
##去掉蛋白长度小于55的基因(55是根据模式植物蒺藜苜蓿5%的蛋白长度阈值设置)
samtools faidx ZM4_V2.pep.fa
awk '$2<55{print $1}' ZM4_V2.pep.fa.fai | grep "t01" |sed 's/\.t01.*//' | sort | uniq > small_protein.txt##转录本对应蛋白信息，只去掉转录本为.t01对应的基因
grep -v -f small_protein.txt ZM4_V2_2.gff3 > ZM4_V2_3.gff3
