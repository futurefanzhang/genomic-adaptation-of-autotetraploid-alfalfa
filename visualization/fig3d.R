setwd("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\kaks_ZM4_V2/")
mydata=read.table("ZM4_V2_Ka_Ks_use.txt",header = F,sep = "\t")
colnames(mydata)=c("quary","hit","blast_pvalue","ka","ks")
mydata=mydata[mydata$ks>0,]
mydata=mydata[mydata$ka>0,]
mydata$kavsks=mydata$ka/mydata$ks  #generate ka/ks
mydata_separated <- mydata %>%
  separate(quary, into = c("gene", "transcript"), sep = ".t")
##use the mean ka/ks of multiple paired gene as the final ka/ks
mydata_result <- mydata_separated %>%
  group_by(gene) %>%
  summarise(mean_kavsks = mean(kavsks), .groups = "drop")
core_gene_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",head = T,sep = "\t")
merged_df <- inner_join(mydata_result, core_gene_info, by = "gene")
#merged_df$sum[is.na(merged_df$sum)] <- "Ungrouped"
merged_df=na.omit(merged_df)
##figure
library(ggplot2)
merged_df$log_kavsks=log(merged_df$mean_kavsks)
merged_df <- merged_df %>%
  mutate(Type = recode(Type, 
                        "one_copy" = "One copy", 
                        "two_copy" = "Two copies", 
                        "three_copy" = "Three copies",
                       "four_copy" = "Four copies"))
                       
merged_df$Type <- factor(merged_df$Type,levels = c("One copy","Two copies","Three copies","Four copies"),ordered = TRUE) 
merged_df$sum <- factor(merged_df$sum,levels = c("Private","Dispensable","Softcore","Core"),ordered = TRUE)
p1=ggplot(data=merged_df,aes(x = sum, y = log_kavsks, fill = Type)) +
  scale_fill_manual(breaks=c("One copy","Two copies","Three copies","Four copies"),
                    values=c("#ABB3CA","#FEC87E","#A2CC00","#FF9477"))+  ##不同类型基因对应的颜色
  geom_boxplot() +
  labs(x = "",y = "log(Ka/Ks)") +
  theme(
    panel.background = element_blank(),
    axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.5, unit = "cm"))
ggsave(p1,filename = "box_core_gene_copy_number_kaks.pdf",width = 7,height = 6,dpi=300)
