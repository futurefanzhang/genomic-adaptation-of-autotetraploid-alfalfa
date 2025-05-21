setwd("E:\\Genome\\Comparative_genomics\\submit_NC\\reply_version1\\lfmm")
climate_gene_info=read.table("240_sample_lfmm_gene.txt",header = T,sep = "\t")
pbs_gene=read.table("E:\\Genome\\Comparative_genomics\\submit_NC\\reply_version1\\PBS\\240sample_PBS_gene.txt",header = T,sep = "\t")
climate_gene_info=rbind(climate_gene_info,pbs_gene)
copy_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",header = T,sep = "\t")
library(dplyr)
climate_gene_info2=copy_info %>% inner_join(climate_gene_info, by=c('gene'='gene'))
result <- climate_gene_info2 %>%
  group_by(Type2) %>%
  summarize(
    count = n(),
    proportion = n() / nrow(climate_gene_info),
    .groups = 'drop'
  )
result

##figure
library(ggplot2)
library(stringr)
climate_fig=read.table("climate_gene_pan_core_copy_number_v2.txt",header = T,sep = "\t")
climate_fig$Gene_group <- str_replace_all(climate_fig$Gene_group, c("one_copy" = "One_copy", "two_copy" = "Two_copy","three_copy" = "Three_copy","four_copy" = "Four_copy","fourfold_core"="Fourfold_core"))

climate_fig$Gene_group=factor(climate_fig$Gene_group,levels = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),ordered = TRUE)
library(ggplot2)
p1 <- ggplot() +
  geom_bar(data=climate_fig,aes(x = Gene_group,y = climate_gene,fill = Gene_group), stat = "identity",  width = 0.8) + 
  scale_fill_manual(
    breaks = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5","#F5D9E6")
  ) +  ## 不同类型基因对应的颜色
  labs(x = "Gene conservation and copy number", y = "Number of climate adaptation related genes") +
  #geom_point(data=climate_fig,aes(x = Copy_number,y = Enrichment*220), color = "red", size = 1) +
  geom_line(data=climate_fig,aes(x = Gene_group,y = Enrichment*400),stat="identity", group=1,color = "red", size = 1) +
  
  scale_y_continuous(sec.axis = sec_axis(~./400, name = "Enrichment")) +
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
ggsave(p1,filename = "bar_climate_gene_pan_core_copy_number_v2.pdf",width = 6,height = 6,dpi=300)
