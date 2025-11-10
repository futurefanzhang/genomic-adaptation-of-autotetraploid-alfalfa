setwd("E:\\Genome\\Comparative_genomics\\submit_NC\\reply_version2\\SV_analysis")
SV_info=read.table("mapped_reads_SV_number.txt",header = T,sep = "\t") #

##mummer4
SV_info1=SV_info[SV_info$group=="mummer4",]
cor(SV_info1$total_sv_count,SV_info1$Mapped.reads)
cor.test(SV_info1$total_sv_count,SV_info1$Mapped.reads)
#minimap2
SV_info2=SV_info[SV_info$group=="minimap2",]
cor(SV_info2$total_sv_count,SV_info2$Mapped.reads)
cor.test(SV_info2$total_sv_count,SV_info2$Mapped.reads)

library(ggplot2)
library(dplyr)
library(ggrepel)
p1=ggplot(data = SV_info1, aes(x = Mapped.reads, y = total_sv_count)) +
  geom_point(aes(color = Species),size = 3) +
  geom_text_repel(aes( color = Species,label = Species),box.padding = unit(0.6, "lines"),point.padding = unit(0.1, "lines"),
                  size = 3, max.overlaps = Inf) +
  geom_smooth(method = lm,formula = y ~ x,se = FALSE,color="black",linewidth = 0.7)+ 
  facet_wrap(~ group, nrow = 1,scales = "free") + 
  annotate(geom="text", x=25000, y=85000, label="Cor=0.99, p=5.13e-16") + 
  labs(title = "",
       x = "Mapped reads",
       y = "Number of SV") +
  scale_color_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D"),
    values = c("#FABB2E","#19B700","#3ea32e","#00868B","#42d4f4","#BF3EFF","#742b99","#A3A500","#dbde43","#6c6e09","#ee0000","#f73e3B","#fa6b6b","#cf4a4a","#940303","#700202","#750e0e","#661a1a"),
    breaks=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17(V5.0)","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2.0_hap2","ZM4_V2.0_hap3","ZM4_V2.0_hap4"))+
  labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "none")
p2=ggplot(data = SV_info2, aes(x = Mapped.reads, y = total_sv_count)) +
  geom_point(aes(color = Species),size = 3) +
  geom_text_repel(aes( color = Species,label = Species),box.padding = unit(0.6, "lines"),point.padding = unit(0.1, "lines"),
                  size = 3, max.overlaps = Inf) +
  geom_smooth(method = lm,formula = y ~ x,se = FALSE,color="black",linewidth = 0.7)+ 
  facet_wrap(~ group, nrow = 1,scales = "free") + 
  annotate(geom="text", x=30000, y=100000, label="Cor=0.95, p=5.31e-10") + 
  labs(title = "",
       x = "Mapped reads",
       y = "Number of SV") +
  scale_color_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D"),
    values = c("#FABB2E","#19B700","#3ea32e","#00868B","#42d4f4","#BF3EFF","#742b99","#A3A500","#dbde43","#6c6e09","#ee0000","#f73e3B","#fa6b6b","#cf4a4a","#940303","#700202","#750e0e","#661a1a"),
    breaks=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17(V5.0)","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2.0_hap2","ZM4_V2.0_hap3","ZM4_V2.0_hap4"))+
  labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "none")

library(cowplot)
figure <- plot_grid(p1, p2,
                    labels = c("a","b"),label_size = 25,ncol = 2, nrow = 1)
ggsave(plot=figure,"mapped_reads_sv_count.pdf", width=15,height=8)
