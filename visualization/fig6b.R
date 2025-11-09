setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\10_USDA_C_N/")
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
merged_fpkm=merged_fpkm %>% inner_join(group_file, by=c('group'='Group'))
genes=read.table("chunxue_MsGDC.txt",header = T,sep = "\t")

wei_genes=merged_fpkm[merged_fpkm$gene %in% genes$gene,]
write.table(wei_genes,"bar_chunxue_gene_exp.txt",row.names = F,col.names = T,sep = "\t",quote = F)

p1= ggplot(data=wei_genes,aes(x = Treatment, y = fpkm, color = gene)) +
  geom_line(aes(linetype = Accession_info,group = interaction( gene, Accession_info)))+facet_grid(. ~ Accession)+
  labs(x = "Treatment",y = "FPKM") +
  theme(
    panel.background = element_blank(),
    axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=14,angle = 45,hjust=1),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.5, unit = "cm"))
ggsave(p1,filename = "bar_chunxue_gene_exp.pdf",width = 10,height = 6,dpi=300)
