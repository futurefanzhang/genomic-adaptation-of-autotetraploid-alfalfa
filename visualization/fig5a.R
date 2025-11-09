mydata=read.table("ZM4_V2_random_70M_gerp_none_zero.txt",header = F,sep="\t") ##原始文件特别大，随机选取7千万行进行统计，由于文件太大，直接在服务器运行R
colnames(mydata)=c("Position","Neutral_rate","RS_score")
gerp_scores <- mydata$RS_score
density_calculation <- density(gerp_scores, bw = "nrd0", n = 512)
density_df <- data.frame(
  x = density_calculation$x,  # x 坐标 (GERP score)
  y = density_calculation$y   # y 坐标 (Density)
)
write.csv(density_df, file = "gerp_density_data.csv", row.names = FALSE)
library(ggplot2)
p1= ggplot(mydata, aes(x = RS_score)) + geom_histogram(aes(y = ..density..),
                                                       colour="black", fill="white") +
  geom_density(alpha = 0.2, fill = "#FF6666") +geom_vline(xintercept = 3, linetype="dashed", color = "blue", size=1.5)+
theme_classic()+labs(x = "GERP score", y = "Density") +
scale_x_continuous(limits = c(-8, 4),breaks = seq(-8,4,1)) +
  theme(
    legend.position = "none",
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(plot=p1,"gerp_distribution.pdf", width=7,height=4)
