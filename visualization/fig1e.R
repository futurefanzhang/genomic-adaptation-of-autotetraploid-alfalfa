mydata=read.table("LTR_all_TE_genome_size.txt",header = T,sep = "\t")
p1=ggplot(mydata, aes(x = genome_size, y = TE_content, color = type2)) +ylim(0,2000)+
  geom_point(size = 3,shape = 1) +  # 绘制点图
  geom_smooth(method = "lm", se = FALSE) +  # 绘制回归线
  labs(title = "",
       x = "Genome size (Mb)",
       y = "TE content (Mb)",
       color = "TE type") +  # 图例标题
  scale_color_manual(values = c("LTR Gypsy & Copia" = "#F9C08A","All TE" = "#A4CB9E")) +  # 自定义颜色
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'),  # 坐标轴标签
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p1,filename = "genome_all_TE_content_figure.pdf",width = 8,height = 5,dpi=300)
