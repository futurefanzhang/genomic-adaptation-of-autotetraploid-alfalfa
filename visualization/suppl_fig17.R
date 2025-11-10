setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\top10_percentage_gene_info/")
library(ggplot2)
library(dplyr)
library(ggrepel)
gene_count=read.table("all_condition_exp_gene.txt",header = T,sep = "\t")  ##top 10% expressed gene and present times
gene_info=read.table("ZM4_V2.bed",header = T,sep = "\t") 
use_info=gene_count %>% inner_join(gene_info, by=c('gene'='gene'))
use_info$SNP=paste(use_info$chromosome,use_info$start, sep="_")
gwas_data=use_info[,c("SNP","Chromosome","start","Count")]
colnames(gwas_data) <- c("ID","CHR", "BP","Count")
gwas_data=gwas_data[order(gwas_data$CHR,gwas_data$BP),]
thres0.01=111.6
##figure
data_cum <- gwas_data %>%
  group_by(CHR) %>%
  summarise(max_bp = as.numeric(max(BP))) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(CHR, bp_add)
gwas_data <- gwas_data %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)
axis_set <- gwas_data %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))
xmin=min(gwas_data$bp_cum)
xmax=max(gwas_data$bp_cum)
##figure
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = Count,
                                  color = as.factor(CHR)),size=1.5) +
  geom_point() +
  geom_hline(yintercept =thres0.01 , color = 'red', linetype = "dashed") +
  scale_x_continuous(label = axis_set$CHR,limits = c(xmin,xmax),expand = expansion(0), breaks = axis_set$center) +
  #scale_y_continuous(expand = expansion(0), limits = c(0, 0.8),breaks = seq(0,0.8,0.2)) +
  ##这块我设置了Fst最大0.8，根据自己的情况改一下
  scale_color_manual(values = rep(c("#F9C08A", "#A4CB9E","#cee2f5","#EDA1A4"), unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome", title = "The frequency of highly expressed genes",
       y = "Count of highly expressed genes") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12,colour = 'black'),
    axis.title.y = element_text(size = 12,colour = 'black'),
    plot.title = element_text(hjust = 0.5,colour = 'black'),
    axis.text.x = element_text(size = 12, angle=30,vjust = 0.5,colour = 'black'),
    axis.text.y = element_text(size = 12,colour = 'black'),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )

##保存png格式
ggsave(plot=manhplot,"count occurence of high expressed genes.png" ,width=16,height=5,units='in',dpi=300)
##保存PDF格式，可以用AI调整细节
ggsave(plot=manhplot,"count occurence of high expressed genes.pdf" ,width=16,height=4)
