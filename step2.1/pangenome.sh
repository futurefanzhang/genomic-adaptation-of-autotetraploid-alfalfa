This script include three major contents, 
###1.PGGB constructed pangenome and calculate sv
###2.SV calculated between two genomes
###3.orthofinder generated pangenes

##1PGGB
perl /data/home/zhangfan/compare_genomoe_medicago/PGGB/change_chr_name.pl change_chr_name.txt Mcalanda.use.fa Mcalanda.pggb.fa ##其中change_chr_name.txt如下
micromamba activate pggb
cat *pggb*fa >in.fa
samtools faidx in.fa
partition-before-pggb -i in.fa -n 19 -o pggb -t 30 -p 90 -s 5k##-n 19是19个单倍型的意思，我这边是11个未分型基因组（不分hap）,2个基因组包括4个hap,所以加起来就是19个, -p 90 核酸片段要有90个碱基一致性， -s 5k 泛基因组
##上一步跑完，就出下面的步骤（来自曹硕的代码），类似：pggb -i in.fa -o output -t 16 -V 'ZM4V2'
pggb -i pggb/in.fa.bf3285f.community.2.fa \
     -o pggb/in.fa.bf3285f.community.2.fa.out \
     -s 5000 -l 25000 -p 90 -c 1 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10M \
     -n 19 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" --resume --threads 30 --poa-threads 30
##gernerate vcf file
vg deconstruct -P ZM4V2 -H '?' -e -a -t 48 in.fa.bf3285f.community.0.fa.out/in.fa.bf3285f.community.0.fa.bf3285f.eb0f3d3.cc60174.smooth.final.gfa > in.fa.bf3285f.community.0.fa.out/0.vcf
perl /data/home/zhangfan/compare_genomoe_medicago/PGGB/pggb/extract_SV_from_vcf.pl 0.vcf chr1_sv.vcf #ref或者alt超过50bp就认为是SV

##2SV, paired genome
conda activate syri_env
nucmer --mum --t 64 -c 1000 --maxgap=500 -p mapping2 /public/home/ac685n40ot/13Medicago/ZM4_V2_hap1_change_chr.fa Mara.chr.fa
delta-filter -1 -i 85 -l 1000 -q mapping2.delta > mapping2.delta-filter
show-coords -THrd mapping2.delta-filter > mapping2.filtered.coords
syri -c mapping2.filtered.coords -d mapping2.delta-filter -r /public/home/ac685n40ot/13Medicago/ZM4_V2_hap1_change_chr.fa -q Mara.chr.fa
perl /public/home/ac685n40ot/13Medicago/Medicago_arabica/filter_syri.pl syri.vcf syri_sv.vcf #过滤syri的结果,保留超过20bp的变异位点
bgzip -c syri_sv.vcf > syri_sv.vcf.gz
tabix -p vcf syri_sv.vcf.gz
minimap2 -a -x asm5 --cs -r2k -t 64 /public/home/ac685n40ot/13Medicago/ZM4_V2_hap1_change_chr.fa Mara.chr.fa  > alignments.sam
samtools sort -m256G -@64 -o alignments.sorted.bam alignments.sam
samtools index alignments.sorted.bam
svim-asm haploid ./ alignments.sorted.bam /public/home/ac685n40ot/13Medicago/ZM4_V2_hap1_change_chr.fa
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz

cuteSV alignments.sorted.bam  /public/home/ac685n40ot/13Medicago/ZM4_V2_hap1_change_chr.fa cuteSV.vcf ./  --threads 12 -s 1 --genotype --report_readid -p -1 -mi 500 -md 500 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
python3 /public/home/ac685n40ot/13Medicago/13genome/diploid_calling.py cuteSV.vcf cuteSV_final.vcf
bgzip -c cuteSV_final.vcf > cuteSV_final.vcf.gz
tabix -p vcf cuteSV_final.vcf.gz
conda activate panpop
bcftools merge -m none -o Merge.vcf variants.vcf.gz cuteSV_final.vcf.gz syri_sv.vcf.gz
sniffles --threads 4 --input alignments.sorted.bam --genotype-vcf Merge.vcf --vcf output_genotypes.vcf #利用sniffles软件把不同软件结果整合在一起
sed -i '/^#/s/Sample\tNULL/Mca_landa/' Mca_landa_genotypes.vcf ##替换sample名称
bgzip -c Mca_landa_genotypes.vcf > Mca_landa_genotypes.vcf.gz
tabix -p vcf Mca_landa_genotypes.vcf.gz
perl /public/home/ac685n40ot/13Medicago/Medicago_arabica/filter_sv_vcf.pl Mara_genotypes.vcf Mara_genotypes_filter.vcf #过滤掉./.的基因型
python /public/home/ac685n40ot/13Medicago/Medicago_arabica/count_frequency.py Mara_genotypes_filter.vcf > SV_count.txt


##3Orthofinder
##version2.5.5
orthofinder -f ./ -oa -M msa -S diamond -t 64 
