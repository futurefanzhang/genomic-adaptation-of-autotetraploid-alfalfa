setwd("E:\\Genome\\PGGB")
diff=read.table("genome_size.txt",header = T,sep = "\t")

diff$Genome=factor(diff$Genome,levels = c("M.truncatula (A17)","ZM4_V2 (hap1)", "ZM4_V2 (Four hap)", "Pan genome"),ordered = TRUE)

library(ggplot2)
p<- ggplot(diff, aes(x=Genome, y=size_gb, fill=Genome)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(
    breaks = c("M.truncatula (A17)","ZM4_V2 (hap1)", "ZM4_V2 (Four hap)", "Pan genome"),
    values = c( "#F9C08A", "#A4CB9E","#cee2f5","#EDA1A4")
  ) +  ## 不同类型基因对应的颜色
  labs(x = "Different genomes", y = "Genome size (Gb)") +
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
ggsave(p,filename = "pangenome_size.pdf",width = 6,height = 6,dpi=300)
