##step0, download the genome used for analysis
##1soybean:ftp://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/glycine_max/dna/
##2Arachis duranensis: https://v1.legumefederation.org/data/v2/Arachis/duranensis/genomes/V14167.gnm2.J7QH/
##3Chickpea:https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=3827&utm_source=assembly
##4Pea:https://urgi.versailles.inra.fr/download/pea/chromosomes/
##5Trifolium pretense: http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/trifolium_pratense/
##6Arabidopsis:ftp://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/arabidopsis_thaliana/dna/
##7Populus trichocarpa:ftp://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/populus_trichocarpa/dna/
##8Vitis vinifera:https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_030704535.1/
###add another 7 genome used in our study
##9 M.arabica, Mca_long, M.lupulina, M.polymorpha,Mru_zhiwusuo, Mtrun_A17(5.0), ZM4_V2

##step1, mkdir ZM4_V2 database, generate mapped maf file 
cd /data/home/zhangfan/compare_genomoe_medicago/GERP_score
for i in $(cat chr_list.txt);do seqkit grep -f $(basename ${i}).txt  ../ZM4_V2/ragtag_chr32_hap1-4.fasta > ZM4_V2_$(basename ${i}).fa; done ##separate chromosome
for i in $(cat chr_list.txt);do lastdb -P4 -uNEAR ZM4_V2_$(basename ${i})_database ZM4_V2_$(basename ${i}).fa; done ##make database for each chromosome
grep ">" Mca.long.fa > chr_info.txt ##提取染色体名称，查看命名方式，然后提取染色体名称到use_chr.txt，去除contig
seqkit grep -f  use_chr.txt Mca.long.fa > Mca.long.use.fa
##seqkit replace -p "(.+?)\s.*" -r '$1' Mpoly.fa > Mpoly.use.fa ##修改ID，仅保留第一个空格之前的ID信息
##seqkit sort -nN Mru.zhiwusuo.use.fa -o Mru.zhiwusuo.use2.fa ##排序染色体
perl ~/compare_genomoe_medicago/ZM4_V2/change_chr_name.pl change_chr_name.txt Mca.long.use.fa Mca.long.use2.fa ##change_chr_name.txt格式如下:Mcalong.chr1, Mcalong.chr2,……,Mcalong.chr8
##change chr name for each genome
##mapping each genome to ZM4_V2 single chromosome
for i in $(cat Species_list.txt);do
last-train -P 20 --revsym --matsym --gapsym -E0.05 -C2 ZM4_V2_Hap1_Chr1_RagTag_database $(basename ${i}) > ZM4_V2_$(basename ${i}).mat
lastal -m50 -E0.05 -C2 -P 20 -p ZM4_V2_$(basename ${i}).mat ZM4_V2_Hap1_Chr1_RagTag_database $(basename ${i}) | last-split -m1 >ZM4_V2_$(basename ${i}).maf
maf-swap ZM4_V2_$(basename ${i}).maf | last-split -m1|maf-swap|maf-sort > ZM4_V2_Hap1_Chr1_RagTag_vs_$(basename ${i}).maf
done
##combine all maf file
sed -i '1 i\##maf version=1 scoring=multiz' ZM4_V2_*_RagTag_vs_*.maf ##添加表头
for a in $(cat Chr_list.txt);do
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara.use.fa.maf ZM4_V2_$(basename ${a})_RagTag_vs_Mca.long.use.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca.maf ZM4_V2_$(basename ${a})_RagTag_vs_Mlup.use.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup.maf ZM4_V2_$(basename ${a})_RagTag_vs_Mpoly.use.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly.maf ZM4_V2_$(basename ${a})_RagTag_vs_Mru.zhiwusuo.use.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru.maf ZM4_V2_$(basename ${a})_RagTag_vs_Mt5.use.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5.maf ZM4_V2_$(basename ${a})_RagTag_vs_aradu.all.chromosome.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu.maf ZM4_V2_$(basename ${a})_RagTag_vs_chickpea.all.chromosome.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea.maf ZM4_V2_$(basename ${a})_RagTag_vs_pea.all.chromosome.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea.maf ZM4_V2_$(basename ${a})_RagTag_vs_soybean.all.chromosome.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean.maf ZM4_V2_$(basename ${a})_RagTag_vs_Trif.all.chromosome.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif.maf ZM4_V2_$(basename ${a})_RagTag_vs_arabidopsis.all.chromosome.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif_arabidopsis.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif_arabidopsis.maf ZM4_V2_$(basename ${a})_RagTag_vs_Populus.all.chromosome.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif_arabidopsis_Populus.maf
multiz ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif_arabidopsis_Populus.maf ZM4_V2_$(basename ${a})_RagTag_vs_vitis.T2T.all.chromosome.use.fa.maf 0 all >ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif_arabidopsis_Populus_vitis.maf
done

##step2, generate phylogenetic tree
##create single copy gene tree, please refer:https://app.yinxiang.com/fx/8da6c1ad-0f59-445c-8fdb-425f1ea183a8
##generate four fold site tree
python filter_MAF_score.py ##过滤maf文件，把ZM4_V2_Hap1_Chr1_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif_arabidopsis_Populus_vitis.maf 文件中不是目标物种Hap1_Chr1_RagTag开头的序列去掉，顺便把score小于20的也去掉
grep "Chr1_1"  /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2_raw.gff3 | grep "CDS" > CDS_Hap1_Chr1.gff3
sed -i 's/Chr1_1/Hap1_Chr1_RagTag/g' CDS_Hap1_Chr1.gff3 ##如果染色体名称不一致的话，替换一下
msa_view filtered.maf --4d --features CDS_Hap1_Chr1.gff3 > 4d-codons.ss
msa_view 4d-codons.ss --in-format SS --out-format SS --tuple-size 1 > 4d-sites.ss
phyloFit --tree all_15_gerp_species.tree.txt --msa-format SS --out-root nonconserved-4d 4d-sites.ss  ##注意，树文件和4d-sites.ss中的物种名称要对应，树结果在nonconserved-4d.mod里面

##step3, calculate GERP value
for a in $(cat Chr_list.txt);do
maf2fasta ZM4_V2_$(basename ${a})_RagTag.fa ZM4_V2_$(basename ${a})_RagTag_vs_Mara_Mca_Mlup_Mpoly_Mru_Mt5_aradu_chickpea_pea_soybean_Trif_arabidopsis_Populus_vitis.maf fasta >ZM4_V2_$(basename ${a})_RagTag_15species.fa
sed "s/Hap1_Chr1/$(basename ${a})/g" 4N_neutral_tree.tree >4N_neutral_tree_$(basename ${a}).tree
/data/home/zhangfan/genome/gerp/gerp/gerpcol -a -f ZM4_V2_$(basename ${a})_RagTag_15species.fa -t 4N_neutral_tree_$(basename ${a}).tree -e $(basename ${a})_RagTag -j  ##这个树用nonconserved-4d.mod里面的结果
done

##step4, add position info
for a in $(cat Chr_list.txt);do
nl ZM4_V2_$(basename ${a})_RagTag_15species.fa.rates > ZM4_V2_$(basename ${a})_RagTag_15species_add_line_number.txt #添加行信息
#source /public/home/wangxu02/anaconda3/bin/activate busco
#python3 /public/home/wangxu02/software/bin/getgaps.py  ZM4_V2.fa > ZM4_V2.gaps.gff3
#grep "Chr1_1" ZM4_V2.gaps.gff3 >ZM4_V2_Hap1_Chr1.gaps.gff3 #提取出Hap1_Chr1的gap
perl /data/home/zhangfan/compare_genomoe_medicago/GERP_score/gerp_add_gap.pl ZM4_V2_$(basename ${a}).gaps.gff3 ZM4_V2_$(basename ${a})_RagTag_15species_add_line_number.txt ZM4_V2_$(basename ${a})_use_gerp.txt
awk -F"\t" '$3>3' ZM4_V2_$(basename ${a})_use_gerp.txt > ZM4_V2_$(basename ${a})_use_gerp_more_than3.txt
done

##step5, calculate gerp value distribution in R
mydata=read.table("ZM4_V2_Hap1_Chr1_use_gerp_none_zero.txt",header = F,sep="\t") ##文件太大，直接在服务器运行R
colnames(mydata)=c("Position","Neutral_rate","RS_score")

p1= ggplot(mydata, aes(x = RS_score)) + geom_histogram(aes(y = ..density..),
                                                       colour="black", fill="white") +
  geom_density(alpha = 0.2, fill = "#FF6666") +geom_vline(xintercept = 3, linetype="dashed", color = "blue", size=1.5)+
theme_classic()+labs(x = "GERP score", y = "Density") +
scale_x_continuous(limits = c(-8, 4),breaks = seq(-8,4,1)) +
  theme(
    legend.position = "none",
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(plot=p1,"gerp_distribution.pdf", width=7,height=4)
mean(mydata$RS_score)
quantile(mydata$RS_score, probs = c(0,0.25,0.27,0.3,0.34,0.37,0.4,0.5,0.75,0.98,1))
