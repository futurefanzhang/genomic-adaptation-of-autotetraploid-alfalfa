setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\DEG_gene_analysis")
file_list <- list.files(pattern = "*.txt") ##put all DEG file in one folder,files similar like ()

# 创建一个空的数据框以存储合并的数据
combined_data <- data.frame()

# 循环读取每个TXT文件并合并
for (file in file_list) {
  # 读取TXT文件
  temp_data <- read.table(file, header = TRUE, sep = "\t")  # 根据需要修改分隔符
  # 合并到主数据框
  combined_data <- rbind(combined_data, temp_data)
}
copy_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",header = T,sep = "\t")
DEG_info=copy_info %>% inner_join(combined_data, by=c('gene'='gene'))
DEG_group <- DEG_info %>%
  group_by(research, stress_type, accession, Type2) %>%
  summarize(
    count = n(),
    .groups = 'drop'
  ) %>%
  group_by(research, stress_type, accession) %>%
  mutate(
    sum=sum(count),
    proportion = count / sum(count)  # 计算比例
  ) 
##filter DEG number greater than 100, because the DEG less than 100 had statistic bias 
DEG_group_filter=DEG_group[DEG_group$sum>=100, ] 
write.table(DEG_group_filter,"all_studies_DEG_proportion_of_gene_group.txt", row.names = F,col.names = T,sep = "\t",quote = F) 
##Check the proportion of fourfold core genes
DEG_group_filter2=DEG_group_filter[DEG_group_filter$Type2=="Fourfold_core",]
##the proportion of fourfold core genes in all stress type of 10_USDA_C_N was different with the rest studies, so we removed it.

##figure
DEG_group_use=DEG_group_filter[DEG_group_filter$research!="10_USDA_C_N",]
##group by stress types
library(stringr)
DEG_group_use$stress_type <- str_replace_all(DEG_group_use$stress_type, c("alkaline_12h" = "Alkaline", "alkaline_24h" = "Alkaline","alkaline_8h" = "Alkaline","Cold_24h" = "Cold","Cold_24h_normal24h"="Cold","Cold_48h"="Cold","Cold_normal24h"="Cold",
                                                                          "no_light" = "Dark", "Salt_and_melatonin" = "Salt and Melatonin","salt" = "Salt","drought" = "Drought","cold" = "Cold"))
##calculate mean and variance
DEG_group_fig=DEG_group_use%>%
group_by(stress_type,Type2) %>%  # 按stress_type分组
  summarize(
    mean_proportion = mean(proportion, na.rm = TRUE),  # 计算均值
    sd_proportion = sd(proportion, na.rm = TRUE),      # 计算标准差
    .groups = 'drop'
  )
DEG_group_fig$sd_proportion[is.na(DEG_group_fig$sd_proportion)] <- 0
library(ggplot2)
DEG_group_fig$Type2 <- str_replace_all(DEG_group_fig$Type2, c("one_copy" = "One_copy", "two_copy" = "Two_copy","three_copy" = "Three_copy","four_copy" = "Four_copy"))

DEG_group_fig$Type2=factor(DEG_group_fig$Type2,levels = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),ordered = TRUE)

p<- ggplot(DEG_group_fig, aes(x=stress_type, y=mean_proportion, fill=Type2)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(
    breaks = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5","#F5D9E6")
  ) +  ## 不同类型基因对应的颜色
  geom_errorbar(aes(ymin=mean_proportion-sd_proportion, ymax=mean_proportion+sd_proportion), width=.2,
                position=position_dodge(.9)) +
  labs(x = "Stress type", y = "DEG gene percentage") +
  theme_classic()+
  theme(
    legend.position="none",
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, colour = 'black', angle = 30,vjust = 0.8),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p,filename = "proportion_deg_stress_type_no_legend_v2.pdf",width = 6,height = 6,dpi=300)
