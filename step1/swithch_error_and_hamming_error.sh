minimap2 -t 40 -ax map-pb --secondary=no /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2_hap4_use.fa muxu.hifi_raw.data.fastq.gz |samtools view -bt 2 /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2_hap4_use.fa.fai - |samtools sort -@ 40 -o ZM4_V2_hap4_hifi.sorted.bam -
samtools index ZM4_V2_hap1_hifi.sorted.bam
samtools sort -@40 -o sorted_ZM4_V2_aln_hifi.bam ZM4_V2_aln_hifi.bam
samtools index sorted_ZM4_V2_aln_hifi.bam

####validate SNP
freebayes -f /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2.fa sorted_ZM4_V2_aln_hifi.bam -r Chr1_1 > freebayes_ZM4_V2_Chr1_1_variants.vcf
##同时处理其他染色体
for i in {2..8}; do
freebayes -f /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2.fa sorted_ZM4_V2_aln_hifi.bam -r Chr${i}_1 > freebayes_ZM4_V2_Chr${i}_1_variants.vcf
done
bcftools filter -e 'QUAL < 30' freebayes_ZM4_V2_Chr1_1_variants.vcf -o freebayes_ZM4_V2_Chr1_1_variants_QUAL30.vcf ##去掉低质量的位点

for i in {2..8}; do
bcftools filter -e 'QUAL < 30' freebayes_ZM4_V2_Chr${i}_1_variants.vcf -o freebayes_ZM4_V2_Chr${i}_1_variants_QUAL30.vcf
whatshap phase --ignore-read-groups -o freebayes_hap1_Chr${i}.phased.vcf --reference=/data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2_hap1.fa --chromosome Chr${i}_1 freebayes_ZM4_V2_Chr${i}_1_variants_QUAL30.vcf sorted_ZM4_V2_aln_ont.bam
done
for i in {1..8}; do
perl /data1/usr/zhangfan/hifiasm/filter_vcf.pl freebayes_hap1_Chr${i}.phased.vcf freebayes_hap1_Chr${i}.phased_GT.vcf  ##仅保留0|1和1|0对应的基因型
done
cat freebayes_hap1_Chr*.phased_GT.vcf|grep -v '#'|grep PS|cut -f1,2,4,5,10 > freebayes_hifi.phase.txt
perl /data1/usr/zhangfan/hifiasm/filter_freebayes_snp.pl freebayes_hifi.phase.txt freebayes_hifi.snp.phase.txt

##reference SNP
nucmer --mum --t 64 -c 1000 --maxgap=500 -p mapping2 /public/home/ac685n40ot/13Medicago/ZM4_V2_hap1_change_chr.fa ZM4_V2.hap2.chr.fa ##Refer:Github:tangerzhang/calc_switchErr
show-snps -Clr mapping2.delta > hap1_hap2.snps
perl /public/home/ac685n40ot/13Medicago/ZM4_V2/hap2/filter_snp_haps.pl hap1_hap2.snps hap1_hap2_filter_chr.snps ##过滤掉不同染色体的比对结果
mv */*_filter_chr.snps ./snp_hap1234/
cat *snps |grep Chr |grep -v home|perl -e 'while(<>){chomp;$_=~s/^\s+//;@t=split(/\s+/,$_);print "$t[13]\t$t[0]\t$t[1]\t$t[13]\t$t[3]\t$t[2]\n" if($t[1] ne "." and $t[2] ne ".")}' > hap1234.hapVar.snps

##calculate switch_error and hamming_error
python calculate_switcherror_hammingerror_v3.py freebayes_hifi.snp.phase.txt hap1234.hapVar.snps
