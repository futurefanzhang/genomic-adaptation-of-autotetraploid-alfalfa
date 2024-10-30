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

#write.table(out,"DESeq2_Cold_normal_11_AllResult.txt", sep = "\t", quote=F,row.names = T,col.names=T)
write.csv(diff,"DESeq2_Cold_normal_11_DEG.csv", row.names = TRUE,quote=F) ##差异基因，差异倍数，P值和FDR值
#write.table(rownames(diff),"DESeq2_Cold_normal_11_DEGlist.txt", sep = "\t",quote=F,row.names=F,col.names=F)

##step3, analysis the fpkm value and generate figures (in R)

setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\2_CAU_cold")
library(DESeq2)
library(dplyr)
countMatrix=read.delim('fpkm.txt',header=T,row.names=1,sep="\t") #FPKM值,整数
countMatrix=round(countMatrix)
group_file=read.delim('group.txt',header=T,sep="\t")  #读取分组文件
###fpkm值
result <-list()
all_fpkm<-list()
groups <- unique(group_file$Group)

# 遍历每个组，选择相应的列并计算均值
for (group in groups) {
  # 获取当前组的列名
  cols <- group_file$Sample[group_file$Group == group]
  # 选择相应的列并计算均值
  col_select <- countMatrix %>% select(all_of(cols))
  mean_values = as.data.frame(apply(col_select, 1, mean))
  # 给结果添加组名
  mean_values$group <- group
  mean_values$gene=rownames(mean_values)
  colnames(mean_values)=c("fpkm","group","gene")
  all_fpkm[[group]]=mean_values
  mean_values=mean_values[mean_values$fpkm>0,]  
  # 确定前5%的阈值
  threshold <- quantile(mean_values$fpkm, 0.95)
  # 筛选出最大值大于或等于阈值的行
  top_5_percent <- mean_values[mean_values$fpkm >= threshold, ]
  result[[group]]=top_5_percent
}
merged_df <- do.call(rbind, result)
merged_fpkm <- do.call(rbind, all_fpkm)
##
hist(mean_values$fpkm,breaks=100,xlim=c(0,1000))
abline(v = quantile(mean_values$fpkm, 0.95), col="red", lwd=3, lty=2)
##
#copy_info=read.table("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\all_gene_pan_core_copy_number_gene_info.txt",head = T,sep = "\t")
#copy_info$Type <- factor(copy_info$Type, levels = c("four_copy","three_copy","two_copy","one_copy"))
#copy_info=copy_info[order(copy_info$gene,copy_info$Type,decreasing = FALSE),]
#copy_info <- copy_info %>%group_by(gene) %>% slice(1)
#write.table(copy_info,"all_gene_pan_core_copy_number_gene_info_expression.txt",row.names = F,col.names = T,sep = "\t",quote = F)
copy_info=read.table("all_gene_pan_core_copy_number_gene_info_expression.txt",header = T,sep = "\t")
result <- copy_info %>%
  group_by( Type2) %>%
  summarize(
    count = n(),
    proportion = n() / nrow(copy_info),
    .groups = 'drop'
  )
result

combine_info=merged_df %>% inner_join(copy_info, by=c('gene'='gene'))
###整体图片
group_file1=group_file %>%group_by(Group) %>% slice(1)
combine_info=combine_info %>% left_join(group_file1, by=c('group'='Group'))
result <- combine_info %>%
  group_by( Type2) %>%
  summarize(
    count = n(),
    proportion = n() / nrow(combine_info),
    .groups = 'drop'
  )
result

##把上面结果整理成表格形式，然后画整体分布图
climate_fig=read.table("fpkm_gene_pan_core_copy_number.txt",header = T,sep = "\t")
climate_fig$Gene_group=factor(climate_fig$Gene_group,levels = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),ordered = TRUE)
p1 <- ggplot() +
  geom_bar(data=climate_fig,aes(x = Gene_group,y = expression_gene,fill = Gene_group), stat = "identity",  width = 0.8) +
  scale_fill_manual(
    breaks = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5","#F5D9E6")
  ) +  ## 不同类型基因对应的颜色
  labs(x = "Gene conservation and copy number", y = "Number of high expressed genes") +
  #geom_point(data=climate_fig,aes(x = Copy_number,y = Enrichment*220), color = "red", size = 1) +
  geom_line(data=climate_fig,aes(x = Gene_group,y = Enrichment*220),stat="identity", group=1,color = "red", size = 1) +
  scale_y_continuous(sec.axis = sec_axis(~./220, name = "Enrichment")) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, angle = 30,vjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.line = element_line(size = 0.5, colour = 'black'), # 轴线
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p1,filename = "bar_fpkm_gene_pan_core_copy_number.pdf",width = 9,height = 6,dpi=300)

library(ggplot2)
library(ggpubr)
combine_info$Type2=factor(combine_info$Type2,levels = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),ordered = TRUE)
p1<-ggplot(combine_info,                         # Draw facet_grid boxplot by split groups
       aes(x = Accession,y = log(fpkm), fill = Accession)) +
  geom_boxplot() + facet_grid(. ~ Type2)+
  labs(x = "Gene conservation and copy number", y = "Log(FPKM)") +
  theme(
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )+ stat_compare_means(comparisons=list(c("Cold sensitive","Cold tolerant")),label = "p.format", label.y = 11)
ggsave(p1,filename = "boxplot_accession_fpkm_gene_pan_core_copy_number.pdf",width = 9,height = 6,dpi=300)
##查看基因在不同条件下的表达量，fpkm值
One_copy_essential=combine_info[combine_info$Type2=="Onefold_essential",]
One_copy_essential$logfpkm=log(One_copy_essential$fpkm)
p1<-ggboxplot(One_copy_essential, x = "Accession", y = "logfpkm",
               color = "Accession", palette = "jco",
               facet.by = "Treatment", short.panel.labs = FALSE)+stat_compare_means(label = "p.format")+
              labs(x = "", y = "Log(FPKM)")
ggsave(p1,filename = "boxplot_accession_fpkm_onefold_essential_treatment.pdf",width = 9,height = 6,dpi=300)
