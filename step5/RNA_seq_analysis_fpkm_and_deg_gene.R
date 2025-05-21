##step1, extract fpkm value (in linux)

for i in $(ls -d */); do echo $i; cd $i; tophat2 --library-type fr-unstranded --read-mismatches 2 -p 20 -G /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2_raw.gff3 -o ./ /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2 *1.fastq.gz *2.fastq.gz; cd ..; done
for i in $(ls -d */); do echo $i; cd $i; samtools view -bh -q 30 -F 4 accepted_hits.bam  |samtools sort >$(basename ${i}).bam; cd ..; done
for i in $(ls -d */); do echo $i; cd $i; htseq-count -r pos -f bam $(basename ${i}).bam /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2_raw.gff3 -t gene -i ID> $(basename ${i})_count_number.txt; cd ..; done
perl /data1/usr/zhangfan/RNA_raw_data/mergeEXP.pl */*count_number.txt all
perl /data1/usr/zhangfan/RNA_raw_data/counts2fpkm.pl all.counts.xls /data/home/zhangfan/compare_genomoe_medicago/ZM4_V2/ZM4_V2_raw.gff3 >fpkm.txt ##保留两位小数

##step2,extract differentially expressed genes (in R)

setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\2_CAU_cold")
library(DESeq2)

countMatrix=read.delim('fpkm.txt',header=T,row.names=1,sep="\t") #FPKM值,整数
countMatrix=round(countMatrix)
group=read.delim('group.txt',header=T,sep="\t")  #读取分组文件
dds <- DESeqDataSetFromMatrix(countMatrix, colData=group, ~ Group)

#过滤掉所有样品中表达量为0的基因
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds,contrast = c("Group","A","B"))   ##填写要比较的分组组合，以B为control，FC=expA/expB
out=cbind(res$log2FoldChange,res$pvalue,res$padj)
rownames(out)=row.names(res)
colnames(out)=c('log2(FoldChange)','Pvalue','FDR')

#===================FDR过滤=============
FCcut=2
FDRcut=0.05
diff=out[(!is.na(res$padj) & res$padj<FDRcut) & abs(res$log2FoldChange)>abs(log2(FCcut)),]
#=======================================

write.csv(diff,"DESeq2_Cold_normal_11_DEG.csv", row.names = TRUE,quote=F) ##差异基因，差异倍数，P值和FDR值
