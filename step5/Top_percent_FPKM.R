setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\23_najing_AU_salt/")
experiment_info="23_najing_AU_salt" ##use one study as example
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

copy_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",header = T,sep = "\t")

##check proportion and export top 10%
for (group in groups) {
combine_info_total=merged_fpkm[merged_fpkm$fpkm>0,]
combine_info_total =combine_info_total %>% inner_join(copy_info, by=c('gene'='gene')) 
combine_info_total=combine_info_total[combine_info_total$group==group,]
threshold <- quantile(combine_info_total$fpkm, 0.9)
top_percent <- combine_info_total[combine_info_total$fpkm >= threshold, ]
top_percent$group2=experiment_info
write.csv(top_percent,paste(experiment_info,"_TOP10%_exp_gene_",group,".csv",sep=""), row.names =F,quote=F)
}

# 筛选出最大值大于或等于阈值的行
result_all=list()
threshold_list=c(0.5,0.75,0.9,0.95,0.99)
a=1
detach("package:dplyr", unload = TRUE)
library(dplyr)
for (group in groups) {
  combine_info_total=merged_fpkm[merged_fpkm$fpkm>0,]
  combine_info_total =combine_info_total %>% inner_join(copy_info, by=c('gene'='gene'))
  combine_info_total=combine_info_total[combine_info_total$group==group,]
  for (i in threshold_list){
    threshold <- quantile(combine_info_total$fpkm, i)
    top_percent <- combine_info_total[combine_info_total$fpkm >= threshold, ]
    result <- top_percent %>%
      group_by(Type2) %>%
      summarize(
        count = n(),
        proportion = n() / nrow(top_percent),
        .groups = 'drop')
    result$total_gene=nrow(top_percent)
    result$group=unique(top_percent$group)
    result$top_percent=i
    result_all[[a]]=result ##小于1不能导出
    a <- a + 1
  }
}
merged_top_info <- do.call(rbind, result_all)
group_file1=group_file %>%group_by(Group) %>% slice(1)
merged_top_info=merged_top_info %>% left_join(group_file1, by=c('group'='Group'))
merged_top_info$group2=experiment_info
write.table(merged_top_info,paste(experiment_info,"All_sample_TOP_percent_gene_exp.txt",sep="_"), sep = "\t", quote=F,row.names = F,col.names=T)
