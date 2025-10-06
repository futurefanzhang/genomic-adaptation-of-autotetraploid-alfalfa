##This script contained two parts of the assembly, ZM4_V1.5 and ZM4_V2.0
##Part 1: ZM4_V1.5

##hifiasm, v0.18.5-r499, generate contig
hifiasm -o muxu_ho -t 100 /public/home/shixiaoya/zf/hifi/muxu.fastq.gz --ul /public/home/shixiaoya/zf/ONT/muxu.ont.fastq.gz --hg-size 800m

##juicer
bwa index /public/home/shixiaoya/zf/Nassembly/hap/juicer/reference/muxu_ho.bp.p_utg.fa
python /public/home/wangxu02/software/juicer/misc/generate_site_positions.py HindIII muxu /public/home/shixiaoya/zf/Nassembly/hap/juicer/reference/muxu_ho.bp.p_utg.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' muxu_HindIII.txt > muxu_HindIII.chrom.sizes       
sh /public/home/wangxu02/software/juicer/CPU/juicer.sh -t 100 -s HindIII -g muxu \
-d /public/home/shixiaoya/zf/Nassembly/hap/juicer \
-D /public/home/wangxu02/software/juicer \
-z /public/home/shixiaoya/zf/Nassembly/hap/juicer/reference/muxu_ho.bp.p_utg.fa \
-p /public/home/shixiaoya/zf/Nassembly/hap/juicer/muxu_HindIII.chrom.sizes \
-y /public/home/shixiaoya/zf/Nassembly/hap/juicer/muxu_HindIII.txt

##3d-dna
/public/home/shixiaoya/soft/3d-dna/run-asm-pipeline.sh -r 0 /public/home/shixiaoya/zf/Nassembly/hap/juicer/reference/muxu_ho.bp.p_utg.fa /public/home/shixiaoya/zf/Nassembly/hap/juicer/aligned/merged_nodups.txt
##manually modification using juicebox
##re3d-dna
/public/home/shixiaoya/soft/3d-dna/run-asm-pipeline-post-review.sh -r /public/home/shixiaoya/zf/Nassembly/hap/re3ddna/muxu_ho.bp.p_utg.final.review.assembly /public/home/shixiaoya/zf/Nassembly/hap/juicer/reference/muxu_ho.bp.p_utg.fa /public/home/shixiaoya/zf/Nassembly/hap/juicer/aligned/merged_nodups.txt

##after these steps, we generated ZM4_V1.5 genome, including the alfalfa haplotype resolved genome with 32 chromosomes (named rename32.fa)

##Part 2: ZM4_V2.0
##after we generate ZM4_V1.5 genome, the genome size is 2.5Gb, it means that part of genome haven't been assembled. so we used ZM4_V1.0 (https://figshare.com/s/fb4ba8e0b871007a9e6c) and ZM4_V1.5 to improve our genome.

##step1: we mapping hifi reads to rename32.fa genome, only keep though high quality reads(MQ>=55)
minimap2 -ax map-hifi /public/agis/zhouyongfeng_group/zhangfan02/genome/new_ref/rename32.fa muxu.hifi.fastq.gz -t 64 > rename32_aln_hifi.sam
samtools view -@ 64 rename32_aln_hifi.sam > mapped_rename32_reads.txt
awk '{if($4>=55) {print $0}}' mapped_rename32_reads.txt | cut -f1,3 > mapped_rename32_high_quality_MQ55_reads.txt #保留质量值超过55的reads,提取reads名称和染色体列，用于后面hifiasm对reads分标签
cut -f1 mapped_rename32_high_quality_MQ55_reads.txt>mapped_rename32_reads_name.txt #提取reads序列的名称

###step2:extract the reads info with MQ<55
samtools view -@ 64 rename32_aln_hifi.sam | cut -f1 >hifi_all_reads_name.txt #提取所有hifi序列名称，包含重复信息
perl /public/agis/zhouyongfeng_group/zhangfan02/genome/hifi_contig/remove_repeat_row.pl hifi_all_reads_name.txt > hifi_unique_reads_name.txt #去掉重复的hifi序列名称，只保留非重复值
perl /public/agis/zhouyongfeng_group/zhangfan02/genome/hifi_contig/compare_file_remove_rep.pl mapped_rename32_reads_name.txt hifi_unique_reads_name.txt unmapped_rename32_reads.txt #提取出来未比对上的hifi的reads名称
seqtk subseq muxu.hifi.fastq.gz unmapped_rename32_reads.txt > unmapped_rename32_raw_hifi_reads.fastq #获取未比对上的fastq文件

##step3: mapping MQ<55 reads to ZM4_V1.0
minimap2 -ax map-hifi zm4_Chr32_long.fa unmapped_rename32_raw_hifi_reads.fastq -t 64 > zm4_aln_hifi.sam
samtools view -@ 64 zm4_aln_hifi.sam > zm4_aln_hifi_flag.txt #去掉标签256,272的序列
awk '{if($4>=55) {print $0}}' zm4_aln_hifi_flag.txt | cut -f1,3 > mapped_zm4_high_quality_MQ55_reads.txt #保留质量值超过55的reads,提取reads名称和染色体列，用于后面hifiasm对reads分标签
perl /public/agis/zhouyongfeng_group/zhangfan02/genome/hifi_contig/name_change1.pl change_chr_info.txt mapped_zm4_high_quality_MQ55_reads.txt zm4_to_rename32_hifi_MQ55_reads.txt ##修改第二个参考基因组和第一个对应关系,change_chr_info.txt两列，第一列旧基因组的染色体，第二列新基因组对应的染色体

##step4: merge MQ>=55 reads in ZM4_V1.5(rename32.fa) and ZM4_V1.0(zm4_Chr32_long.fa)
cat zm4_to_rename32_hifi_MQ55_reads.txt mapped_rename32_high_quality_MQ55_reads.txt > all_hifi_reads_use.txt #合并两个版本的比对结果，然后用于下一步
perl /public/agis/zhouyongfeng_group/zhangfan02/genome/hifi_contig/name_change2.pl change_chr_hap.txt all_hifi_reads_use.txt  mapped_zm4_rename32_MQ55_hap_info.txt ##把32个染色体名称换成hap1-4,统计每个hap对应的reads信息，change_chr_hap.txt包含两列，第一列染色体，第二列hap1-4信息
perl /public/agis/zhouyongfeng_group/zhangfan02/genome/hifi_contig/name_change3.pl change_rename32_chr_to_group.txt all_hifi_reads_use.txt mapped_zm4_rename32_MQ55_hifi_use.txt #change_rename32_chr_to_group.txt包含两列，第一列是染色体，第二列是连锁群

##re-run hifiasm, 0.19.7-r598
#pwd, /data1/usr/zhangfan/hifiasm/
hifiasm -o muxu_four_hap_rename32_new_group -t 80 -D 10 --ul muxu.ont.longerthan50kb.fastq muxu.hifi_raw.data.fastq.gz -5  mapped_zm4_rename32_MQ55_hifi_use.txt --n-hap 4 --hom-cov 120
awk '/^S/{print ">"$2;print $3}' muxu_four_hap_rename32_new_group.bp.hap1.p_ctg.gfa > muxu_four_hap_rename32_new_group.bp.hap1.p_ctg.fa #same for rest haps
seqkit stats muxu_four_hap_rename32_new_group.bp.hap*.p_ctg.fa -a

##using ragtag to generate ZM4_V2.0 genome, single8.fa was the hap1 genome of ZM4_V1.5
python /public/home/wangxu02/anaconda3/envs/ragtag/bin/ragtag.py scaffold /public/agis/zhouyongfeng_group/zhangfan02/genome/new_ref/single8.fa /public/agis/zhouyongfeng_group/zhangfan02/genome/hifi_contig/muxu_four_hap_rename32_new_group.bp.hap1.p_ctg.fa -o hap1.scaffold.fa ##same for rest haps
