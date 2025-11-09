data=read.table("GO_BP_figure_v2.txt",header = T,sep="\t")

df <- data[order(data$Fold.Enrichment),]
library(stringr)
library(ggplot2)
p1<-ggplot(df,aes(x = Fold.Enrichment, y = bioprocess)) +        # x 轴用GeneRatio, y轴用GO或KEGG注释
  geom_point(aes(color = -log10(PValue), size = Count), pch = 19) +            # 颜色用p.adjust，size设为用count
  scale_size_continuous( range = c(3,7))+
  xlim(1,3)+labs(title = "", x = "Fold enrichment",y = "")+
  scale_color_steps(low = "orange", high = "red") +  
  #theme_classic() +
  # theme(axis.text.x = element_text(angle = 45,hjust=0.5)) +
  theme_bw() +
  theme(#panel.grid=element_blank(), #去掉网格线
    axis.title.x = element_text(color="black",size = 12),
    axis.title.y = element_text(color="black",size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black",size = 12, vjust = 0.5),
    axis.text.y = element_text(color="black",size = 12))

ggsave(plot=p1,"GO_enrichment_four_hap_deleterious_v2.pdf" ,width=8,height=5)
