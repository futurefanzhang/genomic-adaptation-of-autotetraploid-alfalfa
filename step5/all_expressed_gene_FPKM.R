setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\8_SASK_salt_rhizobium/") ##one example
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
merged_fpkm=merged_fpkm[merged_fpkm$fpkm>0,] ##only use fpkm value larger than 0
merged_fpkm=merged_fpkm %>% inner_join(group_file, by=c('group'='Group'))
mean_fpkm <- merged_fpkm %>%
  group_by(Treat_group, gene) %>%
  summarise(mean_fpkm = mean(fpkm, na.rm = TRUE))
mean_fpkm$group="Abiotic stress"

write.table(mean_fpkm,"all_gene_expression_FPKM_value_8_SASK_salt_rhizobium.txt",row.names = F,col.names = T,sep = "\t",quote = F)
