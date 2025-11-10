setwd("E:\\Genome\\Comparative_genomics\\submit_NC\\reply_version2\\reads_mapping_info/")
df=read.table("ONT_hifi_reads.txt",header = T,sep = "\t")
library(dplyr)

df_processed <- df 
##绘图
library(ggprism)
library(ggbreak)
p <- ggplot(data = df_processed, aes(x = Haps, y = Count, fill = Type)) +
  geom_col() +
  facet_grid(. ~ Sequencing) +
  # 将 expand 移到这里，因为它对两个Y轴都适用
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Haplotypes", y = "Reads count") +
  
  # --- 主题设置 ---
  theme_prism(
    palette = "flames",
    base_fontface = "plain",
    base_family = "serif",
    base_size = 16,
    base_line_size = 0.8
  ) +
  
  theme(
    # 将图例移到右侧
    legend.position = "right",
    text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),   # 坐标轴标题 (X和Y)
    axis.text = element_text(colour = "black"),    # 坐标轴刻度文字 (X和Y)
    strip.text = element_text(colour = "black"),   
    # 1. 彻底移除右侧Y轴的所有元素
    #    ggbreak 会创建一个 secondary Y axis，我们需要把它关掉
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(), # <-- 关键补充：移除右侧轴线
    
    # 3. 将坐标轴线和刻度线颜色设为黑色
    axis.line = element_line(colour = "black"), # 同时设置X和Y轴线
    axis.ticks = element_line(colour = "black") # 同时设置X和Y轴刻度线
  )
p2 <- p +
  scale_y_break(c(31000, 1300000), 
                scales = "free", # 让上下两部分的刻度更合理
                ticklabels = c(1300000, 1400000,1500000)) # 可选：手动指定上部刻度
ggsave(plot=p2,"reads_mapped_count.pdf", width=8,height=5, device = cairo_pdf )
