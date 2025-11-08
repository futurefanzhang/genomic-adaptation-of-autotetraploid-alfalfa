setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\DEG_gene_analysis")
file_list <- list.files(pattern = "*DEG.txt")

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
write.table(DEG_info,"all_DEG_FC_six_group.txt", row.names = F,col.names = T,sep = "\t",quote = F) 

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

##figure
DEG_group_filter=read.table("E:\\Genome\\Comparative_genomics\\submit_NC\\reply_version3\\supply_table\\Supplementary Data\\supply_table7_Information on differentially expressed genes among different research groups.txt",header = T,sep = "\t")
DEG_group_use=DEG_group_filter[DEG_group_filter$research!="7_USDA_C_N",] ##the nitrogen and CO2 levels can't be defined as stress, so we removed it.
##group by stress types
library(stringr)
DEG_group_use$stress_conditions <- str_replace_all(DEG_group_use$stress_conditions, c("alkaline_12h" = "Alkaline", "alkaline_24h" = "Alkaline","alkaline_8h" = "Alkaline","Cold_24h" = "Cold","Cold_normal24h"="Cold","Cold_48h"="Cold","Cold_24h_normal24h"="Cold",
                                                                          "no_light" = "Dark", "Salt_and_melatonin" = "Salt and Melatonin","salt" = "Salt","drought" = "Drought","cold" = "Cold"))
DEG_group_use$gene_groups <- str_replace_all(DEG_group_use$gene_groups, c("one_copy" = "One_copy", "two_copy" = "Two_copy","three_copy" = "Three_copy","four_copy" = "Four_copy"))

DEG_group_use$gene_groups=factor(DEG_group_use$gene_groups,levels = c("Unique_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Tetra_copy_core"),ordered = TRUE)

##calculate mean and variance
DEG_group_fig=DEG_group_use%>%
group_by(stress_conditions,gene_groups) %>%  # 按stress_type分组
  summarize(
    mean_proportion = mean(proportion, na.rm = TRUE),  # 计算均值
    sd_proportion = sd(proportion, na.rm = TRUE),      # 计算标准差
    .groups = 'drop'
  )
DEG_group_fig$sd_proportion[is.na(DEG_group_fig$sd_proportion)] <- 0
library(ggplot2)

##box plot
p_boxplot_points <- ggplot(DEG_group_use, 
                           # 在全局 aes() 中同时定义 fill 和 color
                           aes(x=stress_conditions, y=proportion, fill=gene_groups, color=gene_groups)) +
  
  # 1. 绘制箱线图 (基本不变)
  geom_boxplot(
    position = position_dodge(0.8),
    width = 0.7,
    color = "black", 
    outlier.shape = NA,
  ) +
  
 
  geom_point(
    # 使用 position_dodge 让点和箱线图对齐
    position = position_dodge(0.8),
    size = 1,           # 点的大小
    shape = 21,         # 使用带边框的圆形 (关键！)
    color = "black",    # 将点的边框设为黑色，更清晰
    stroke = 0.5        # 边框的粗细
  ) +

  # 只需要一个 scale_fill_manual，ggplot会自动应用于 color
  scale_fill_manual(
    name = "Gene Groups",
    breaks = c("Unique_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Tetra_copy_core"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5","#F5D9E6")
  ) +
  # 再加一个 scale_color_manual 确保点的颜色也一样
  scale_color_manual(
    name = "Gene Groups",
    breaks = c("Unique_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Tetra_copy_core"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5","#F5D9E6")
  ) +
  
  labs(x = "Stress Conditions", y = "DEG gene percentage") +
  theme_classic() +
  theme(
    legend.position="none",
    axis.title = element_text(size = 14, colour = 'black'),
    axis.text = element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm")
  )


ggsave(p_boxplot_points,filename = "proportion_deg_stress_type_no_legend_v3.pdf",width = 6,height = 6,dpi=300)
