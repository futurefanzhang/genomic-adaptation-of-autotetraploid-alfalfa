setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\all_expressed_gene_FPKM_value")
library(dplyr)
library(purrr)  # 
file_list <- list.files(path = "./", pattern = "^all_gene_expression", full.names = TRUE)
# 批量读取所有 .txt 文件，将它们存储在一个列表中，并选择所需的列
data_list <- lapply(file_list, function(file) {
  read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    select(Treat_group, gene, mean_fpkm,group)
})
merged_exp_info <- bind_rows(data_list)
#merged_exp_info=read.table("all_gene_expression_FPKM_value_Long_new_salt.txt",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
##Abiotic stress
merged_exp_info_ab=merged_exp_info[merged_exp_info$group=="Abiotic stress",]
merged_exp_info_ab$Treat_group[merged_exp_info_ab$Treat_group == "Background"] <- "Normal"
merged_exp_info_ab$Treat_group[merged_exp_info_ab$Treat_group == "Cold"] <- "Stress"
mean_fpkm_ab <- merged_exp_info_ab %>%
  group_by(Treat_group, gene) %>%
  summarise(mean_fpkm = mean(mean_fpkm, na.rm = TRUE))
mean_fpkm_ab$Treat_group[mean_fpkm_ab$Treat_group == "Stress"] <- "Abiotic stress"
mean_fpkm_ab$Treat_group[mean_fpkm_ab$Treat_group == "Stree"] <- "Abiotic stress"
##biotic stress
merged_exp_info_bi=merged_exp_info[merged_exp_info$group=="Biotic stress",]
mean_fpkm_bi <- merged_exp_info_bi %>%
  group_by(Treat_group, gene) %>%
  summarise(mean_fpkm = mean(mean_fpkm, na.rm = TRUE))
mean_fpkm_bi$Treat_group[mean_fpkm_bi$Treat_group == "Stress"] <- "Biotic stress"
mean_fpkm_bi$Treat_group[mean_fpkm_bi$Treat_group == "Stree"] <- "Biotic stress"
mean_fpkm_all=rbind(mean_fpkm_ab,mean_fpkm_bi)
mean_fpkm_all <- mean_fpkm_all %>%
  group_by(Treat_group, gene) %>%
  summarise(mean_fpkm = mean(mean_fpkm, na.rm = TRUE))
combine_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",head = T,sep = "\t")
mean_fpkm_all=mean_fpkm_all %>% inner_join(combine_info, by=c('gene'='gene'))
##figure
mean_fpkm_fig=mean_fpkm_all
mean_fpkm_fig$Type2 <- str_replace_all(mean_fpkm_fig$Type2, c("Onefold_essential"="Unique_essential","one_copy" = "One_copy", "two_copy" = "Two_copy","three_copy" = "Three_copy","four_copy" = "Four_copy","Fourfold_core"="Tetra_copy_core"))

mean_fpkm_fig$Type2=factor(mean_fpkm_fig$Type2,levels = c("Unique_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Tetra_copy_core"),ordered = TRUE)
mean_fpkm_fig$log_fpkm=log2(mean_fpkm_fig$mean_fpkm+1)
library(stringr)
write.table(mean_fpkm_fig,"all_exp_gene_fpkm.txt",row.names = F,col.names = T,sep = "\t",quote = F)
p<- ggplot(mean_fpkm_fig, aes(x=Type2, y=log_fpkm, color=Type2)) + 
  geom_boxplot() +
  scale_color_manual(
    breaks = c("Unique_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Tetra_copy_core"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5","#F5D9E6")
  ) +  ## 不同类型基因对应的颜色
  labs(x = "Groups", y = "Log2(FPKM+1)") +
  theme_classic()+facet_grid(. ~ Treat_group)+
  theme(
    legend.position = "right",
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p,filename = "all_exp_gene_fpkm.pdf",width = 7,height = 6,dpi=300)

##Which genes high exp among tetra-copy core genes
tetra_copy_core=mean_fpkm_all[mean_fpkm_all$Type2=="Fourfold_core",]
mean_fpkm_4core <- tetra_copy_core %>%
  group_by(gene) %>%
  summarise(mean_fpkm = mean(mean_fpkm, na.rm = TRUE))
four_pair=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\combine_four_hap_gene1234_paired_gene.txt",header = T,sep = "\t")
exp_four_pair=four_pair %>% left_join(mean_fpkm_4core, by=c('gene_name'='gene'))
exp_four_pair=exp_four_pair[,c("gene_name","pair","mean_fpkm")]
filtered_exp_four_pair <- exp_four_pair %>%  ##去掉不表达的fourfold_core gene的信息
  group_by(pair) %>%
  filter(!all(is.na(mean_fpkm)))
filtered_exp_four_pair <- filtered_exp_four_pair %>% mutate(mean_fpkm = ifelse(is.na(mean_fpkm), 0, mean_fpkm))
filtered_exp_four_pair <- filtered_exp_four_pair %>%
  group_by(pair) %>%
  arrange(desc(mean_fpkm), .by_group = TRUE) %>%
  mutate(rank = row_number()) ##order four haps by the FPKM value 
##figure
filtered_exp_four_pair$logfpkm=log2(filtered_exp_four_pair$mean_fpkm+1)
grouped_means <- filtered_exp_four_pair %>%
  group_by(rank) %>%
  summarise(mean_fpkm_mean = mean(mean_fpkm))
grouped_means 
pairwise_results <- pairwise.t.test(filtered_exp_four_pair$mean_fpkm, filtered_exp_four_pair$rank, p.adjust.method = "bonferroni")
pairwise_results
filtered_exp_four_pair$rank<-as.factor(filtered_exp_four_pair$rank)

df.summary <- filtered_exp_four_pair %>%
  group_by(rank) %>%
  summarise(
    sd = sd(logfpkm, na.rm = TRUE),
    len = mean(logfpkm)
  )
df.summary
##mean and sd figure
pairwise_results <- pairwise.t.test(filtered_exp_four_pair$logfpkm, filtered_exp_four_pair$rank, p.adjust.method = "bonferroni")
pairwise_results
df.summary$significance="p < 2e-16"
p <- ggplot(
  df.summary, 
  aes(x = rank, y = len, ymin = len-sd, ymax = len+sd)
)+ geom_errorbar(width = 0.3,color="#F9C08A") + geom_point(size = 2) +
  geom_segment(aes(x = 0.9, xend = 4, y = 4)) +  # 添加横线
  geom_text(x = 1.25, y = 4.1, label = "p < 2e-16", size = 4) +  # 在横线上添加星号
  geom_segment(aes(x = 0.9, xend = 3, y = 3.8)) +  # 添加横线
  geom_text(x = 1.25, y = 3.9, label = "p < 2e-16", size = 4) +  # 在横线上添加星号
  geom_segment(aes(x = 0.9, xend = 2, y = 3.6)) +  # 添加横线
  geom_text(x = 1.25, y = 3.7, label = "p < 2e-16", size = 4) +  # 添加P值
  labs(title = "",
       x = "The expression levels ranked in descending order",
       y = "Log2(FPKM+1)") +theme_classic()+
  theme(
    legend.position = "none",
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p,filename = "all_tetra-copy_core_gene_mean_fpkm_order1234_V2.pdf",width = 5,height = 6,dpi=300)
