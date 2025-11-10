setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\DEG_fold_change")
diff=read.table("all_DEG_FC_six_group.txt",header = T,sep = "\t")
all_studied_research=read.table("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\DEG_gene_analysis\\all_studies_DEG_proportion_of_gene_group.txt",header = T,sep = "\t")
diff$Type2 <- str_replace_all(diff$Type2, c("one_copy" = "One_copy", "two_copy" = "Two_copy","three_copy" = "Three_copy","four_copy" = "Four_copy"))

diff$Type2=factor(diff$Type2,levels = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),ordered = TRUE)
diff=diff[diff$research %in% all_studied_research$research,]
p1<-ggplot(diff,                         # Draw facet_grid boxplot by split groups
           aes(x = Type2,y = log(abs(log2.FoldChange.)), fill = Type2)) +
  geom_boxplot() +
  scale_fill_manual(
    breaks = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5","#F5D9E6")
  ) +  ## 不同类型基因对应的颜色
  labs(x = "Gene conservation and copy number", y = "Log(|Log2(FoldChange)|)") +
  theme_classic()+
  theme(
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14),
    legend.position = "none" ,
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p1,filename = "boxplot_diff_exp_all_type.pdf",width = 15,height = 6,dpi=300)
