##suppl fig20
df=read.table("suppl_fig20.txt",header = T,sep = "\t")
df$lines=factor(df$lines,levels = c("WT", "MsG1", "MsG2", "MsGS1","MsGS10"),ordered = TRUE)

df1=df[df$phenotype=="Fresh weight",]
p1= ggplot(data = df1, aes(x = lines, y = value,fill = lines)) +
  
  # 绘制均值柱状图
  # stat_summary() 会自动计算均值并绘制几何图形（这里是"bar"）
  stat_summary(fun = "mean", 
               geom = "bar", 
               
               color = "black",    # 柱子边框为黑色
               width = 0.6) +     # 柱子宽度
  
  # 叠加原始数据点
  
  # 按照 "chlorophy" 列进行分面
  geom_point(size = 3.5,       # 点的大小，可以适当调整
             shape = 1,         # shape = 1 代表空心圆圈
             colour = "black",  # 点的颜色
             stroke = 1.1) +    # 点的线条粗细，让空心圆更明显
  
  scale_fill_manual(
    breaks = c("WT", "MsG1", "MsG2", "MsGS1","MsGS10"),
    values = c("gray50", "#74A9CF","#034E7B", "#F7A8B8","#C51B7D")
  ) +  ## 不同类型基因对应的颜色
  # 添加坐标轴标签
  labs(x = "Lines", y = "Fresh weight (g)") +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14, colour = 'black'),
    axis.text = element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"),
    # 额外美化一下分面的标题，使其更清晰
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
  )
df2=df[df$phenotype=="Dry weight",]
p2= ggplot(data = df2, aes(x = lines, y = value,fill = lines)) +
  
  # 绘制均值柱状图
  # stat_summary() 会自动计算均值并绘制几何图形（这里是"bar"）
  stat_summary(fun = "mean", 
               geom = "bar", 
               
               color = "black",    # 柱子边框为黑色
               width = 0.6) +     # 柱子宽度
  
  # 叠加原始数据点
  
  
  # 按照 "chlorophy" 列进行分面
  geom_point(size = 3.5,       # 点的大小，可以适当调整
             shape = 1,         # shape = 1 代表空心圆圈
             colour = "black",  # 点的颜色
             stroke = 1.1) +    # 点的线条粗细，让空心圆更明显
  
  scale_fill_manual(
    breaks = c("WT", "MsG1", "MsG2", "MsGS1","MsGS10"),
    values = c("gray50", "#74A9CF","#034E7B", "#F7A8B8","#C51B7D")
  ) +  ## 不同类型基因对应的颜色
  # 添加坐标轴标签
  labs(x = "Lines", y = "Dry weight (g)") +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14, colour = 'black'),
    axis.text = element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"),
    # 额外美化一下分面的标题，使其更清晰
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
  )
df3=df[df$phenotype=="Plant height",]
p3= ggplot(data = df3, aes(x = lines, y = value,fill = lines)) +
  
  # 绘制均值柱状图
  # stat_summary() 会自动计算均值并绘制几何图形（这里是"bar"）
  stat_summary(fun = "mean", 
               geom = "bar", 
               
               color = "black",    # 柱子边框为黑色
               width = 0.6) +     # 柱子宽度
  
  # 叠加原始数据点
  
  # 按照 "chlorophy" 列进行分面
  geom_point(size = 3.5,       # 点的大小，可以适当调整
             shape = 1,         # shape = 1 代表空心圆圈
             colour = "black",  # 点的颜色
             stroke = 1.1) +    # 点的线条粗细，让空心圆更明显
  
  scale_fill_manual(
    breaks = c("WT", "MsG1", "MsG2", "MsGS1","MsGS10"),
    values = c("gray50", "#74A9CF","#034E7B", "#F7A8B8","#C51B7D")
  ) +  ## 不同类型基因对应的颜色
  # 添加坐标轴标签
  labs(x = "Lines", y = "Plant height (cm)") +
  
  # 应用您指定的绘图偏好
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14, colour = 'black'),
    axis.text = element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"),
    # 额外美化一下分面的标题，使其更清晰
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
  )
library(cowplot)
#combine
figure <- plot_grid(p1, p2,p3,
                    labels = c("A","B","C"),label_size = 25,ncol = 3, nrow = 1)

ggsave(figure,filename = "suppl_fig20.pdf",width = 10,height = 6,dpi=300)

##suppl fig21
mydata=read.table("suppl_fig21.txt",header = T,sep = "\t")
library(ggplot2)
mydata$Line=factor(mydata$Line,levels = c("WT", "MsG1", "MsG2", "MsGS1","MsGS10"),ordered = TRUE)
mydata1 <- mydata %>%
  group_by(tissue, Line, environment) %>%
  summarise(
    expression_mean = mean(expression),
    sd = sd(expression),
    .groups = 'drop' # 推荐加上.groups='drop'
  )
p <- ggplot() + 
  
  # --- 第1层：使用汇总数据 mydata 绘制均值柱状图 ---
  geom_bar(data = mydata1, 
           aes(x = tissue, y = expression_mean, fill = Line),
           stat = "identity", 
           color = "black", 
           position = position_dodge(width = 0.9)) +
  
  geom_point(data = mydata,
             aes(x = tissue, y = expression, group = Line), # 必须指定 group = Line
             position = position_dodge(width = 0.9),      # ★ 关键：让点在组间躲避
             size = 2.5,       # 点的大小
             shape = 1,        # shape = 1 代表空心圆圈
             colour = "black", # 点的颜色
             stroke = 1) +     # 点的线条粗细
  
  
  scale_fill_manual(
    breaks = c("WT", "MsG1", "MsG2", "MsGS1","MsGS10"),
    values = c("gray50", "#74A9CF","#034E7B", "#F7A8B8","#C51B7D")
  ) +
  labs(x = "Tissue", y = "Relative expression") +
  theme_classic() +
  facet_grid(. ~ environment) +
  theme(
    legend.position="right",
    axis.title = element_text(size = 14, colour = 'black'),
    axis.text = element_text(size = 14, colour = 'black'),
    axis.text.x = element_text(size = 14, colour = 'black', angle = 30, vjust = 0.8),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"),
    # 可以给分面加个边框，让区分更清晰
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )
ggsave(p,filename = "supply_fig21.pdf",width = 10,height = 6,dpi=300)


##suppl fig22
df=read.table("suppl_fig22.txt",header = T,sep = "\t")
df$lines=factor(df$lines,levels = c("WT", "MtG1", "MtG5", "MtGS11","MtGS12"),ordered = TRUE)

p1= ggplot(data = df, aes(x = lines, y = value,fill = lines)) +
  
  # 绘制均值柱状图
  # stat_summary() 会自动计算均值并绘制几何图形（这里是"bar"）
  stat_summary(fun = "mean", 
               geom = "bar", 
               
               color = "black",    # 柱子边框为黑色
               width = 0.6) +     # 柱子宽度
  
  # 叠加原始数据点
  
  
  # 按照 "chlorophy" 列进行分面
  geom_point(size = 3.5,       # 点的大小，可以适当调整
             shape = 1,         # shape = 1 代表空心圆圈
             colour = "black",  # 点的颜色
             stroke = 1.1) +    # 点的线条粗细，让空心圆更明显
  facet_wrap(~ df$chlorophyll.type) +
  scale_fill_manual(
    breaks = c("WT", "MtG1", "MtG5", "MtGS11","MtGS12"),
    values = c("gray50", "#74A9CF","#034E7B", "#F7A8B8","#C51B7D")
  ) +  ## 不同类型基因对应的颜色
  # 添加坐标轴标签
  labs(x = "Lines", y = "Chlorophyll Content (mg/g)") +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14, colour = 'black'),
    axis.text = element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"),
    # 额外美化一下分面的标题，使其更清晰
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )
ggsave(p1,filename = "suppl_fig22.pdf",width = 10,height = 6,dpi=300)
