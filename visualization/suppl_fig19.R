setwd("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info")
gene_info=read.table("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\top10_percentage_gene_info\\ZM4_V2.bed",header = T,sep = "\t") 
copy_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",header = T,sep = "\t")
library(dplyr)
combine=gene_info %>% inner_join(copy_info, by=c('gene'='gene'))
##fourfold core genes
combine1=combine[combine$Type2=="Fourfold_core",]
write.table(combine1,"ZM4_gene_Fourfold_core_gene.txt", row.names = F,quote=F,sep = "\t")
##each type gene density
library(dplyr)
library(tidyr)
# 定义区间大小
interval_size <- 1000000  # 1Mb

# 创建一个函数来计算区间
compute_interval <- function(start, interval_size) {
  interval_start <- floor(start / interval_size) * interval_size + 1
  interval_end <- interval_start + interval_size - 1
  return(list(interval_start = interval_start, interval_end = interval_end))
}

# 应用函数生成区间列
use_info <- combine1 %>%
  rowwise() %>%
  mutate(interval_info = list(compute_interval(start.x, interval_size))) %>%
  unnest_wider(interval_info)

# 统计每个染色体和区间内基因类型的出现次数
result_fourfold_core <- use_info %>%
  group_by(Chromosome.x, interval_start, interval_end) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(Chromosome.x, interval_start, interval_end)

write.table(result_fourfold_core,"ZM4_gene_density_Fourfold_core_gene.txt", row.names = F,quote=F,sep = "\t")

# 应用函数生成区间列
use_info <- combine %>%
  rowwise() %>%
  mutate(interval_info = list(compute_interval(start.x, interval_size))) %>%
  unnest_wider(interval_info)

# 统计每个染色体和区间内基因类型的出现次数
result_all_gene <- use_info %>%
  group_by(Chromosome.x, interval_start, interval_end) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(Chromosome.x, interval_start, interval_end)
write.table(result_all_gene,"ZM4_gene_density_all_gene.txt", row.names = F,quote=F,sep = "\t")
##fourfold core enriched (fourfold/all_gene)
combine_gene_density=result_all_gene %>% inner_join(result_fourfold_core, by=c('Chromosome.x'='Chromosome.x','interval_start'='interval_start'))
combine_gene_density$fourcopy_core_density=(combine_gene_density$count.y)/sum(combine_gene_density$count.y)
combine_gene_density$all_gene_density=(combine_gene_density$count.x)/sum(combine_gene_density$count.x)
combine_gene_density$fourcopyc_core_enriched=combine_gene_density$fourcopy_core_density/combine_gene_density$all_gene_density
##heatmap
chr_info=read.table("ZM4_V2_chr_length.txt",header = T,sep = "\t")
mydata2=combine_gene_density
mydata2$fourcopyc_core_enriched[mydata2$fourcopyc_core_enriched > 2] <- 2
colnames(chr_info)=c("chr_name","chr_length")
####绘图####R版本4.1.2，否则会有边缘重影
library(ggplot2)
write.table(mydata2,"fourfold_core_gene_density_enriched.txt",row.names = F,col.names = T,sep = "\t",quote = F)

p1=ggplot(data=chr_info,aes(x=chr_name,y=chr_length/1000000)) + #载入数据
  ####绘制染色体竖向长条####
geom_col(position=position_nudge(x = 0), width = 0.5, color = "gray60",fill="white") + 
  #绘制染色体长度的长条，往左侧偏移0，宽度为0.1,黑色
  scale_y_continuous(trans = "reverse")+ #把纵坐标进行反向排列（0在最上面）
  ####绘制SNP位置小横线####
geom_tile(data = mydata2, aes(x = factor(Chromosome.x), y = interval_start/1000000, fill = fourcopyc_core_enriched), color = NA,
          size = 0.5, width = .5,position=position_nudge(x = 0))+ 
  #绘制snp横线，每个snp一条横线
  scale_fill_gradient(low = "#440154", high = "#FDE725",limits = c(0, 2), name  = "Tetra-copy core genes enrichment") +#修改不同snp类型的颜色
  theme(plot.margin = unit(x = c(3, 3, 3, 3), units = "mm"), #图像边缘留空白
        panel.background = element_blank(),text = element_text(),legend.position = c(0.3,0.1),legend.direction = "horizontal") + #legend.position=c(0.91,0.95), 0.91 represent left and right, 0.95 represent up and down
  labs(x = 'Chromosome', y = 'Physical Position (Mb)',
       title = '')+theme(axis.line = element_line(colour = "black", size = 0.5),
                         legend.title = element_text(color = "black", size = 12), 
                         legend.text = element_text(color = "black", size = 12),
                         axis.text.y = element_text(size = 12, color = "black"),
                         axis.text.x = element_text(size = 12, color = "black",angle = 30,vjust = 0.5),
                         axis.title.x=element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size=12),
                         axis.title.y=element_text(size=12))

ggsave(plot = p1,file = "fourfold_core_gene_density_enriched_v2.pdf",width = 16,height=6) 
