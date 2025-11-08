setwd("E:\\Genome\\genome_SV")
SV_info=read.table("SV_count.txt",header = T,sep = "\t") #
SV_use=na.omit(SV_info) ##去掉非SV的信息
SV_result <- SV_use %>%
  group_by(Species, SV_type) %>%  # 按物种和SV类型分组
  summarise(total_count = sum(Count, na.rm = TRUE)) %>%  # 对count进行求和
  ungroup()  # 解除分组
SV_summary=SV_result  %>% group_by(Species) %>% summarise(total_count = sum(total_count, na.rm = TRUE))
write.table(SV_summary,"SV_summary_18haps.txt",row.names = F,col.names = T,sep = "\t",quote = F)
##figure
SV_result$SV_type <- factor(SV_result$SV_type,levels = c("TRANS","INV","CNV","DUP","DEL","INS"),ordered = TRUE) 
SV_result$Species <- factor(SV_result$Species,levels = c("Mara","Mlup","Mpoly","Mru_landa","Mru_zhiwusuo","Mt_HM078","Mt_A17","Mt_R108","Mca_landa","Mca_long","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2_hap2","ZM4_V2_hap3","ZM4_V2_hap4"),ordered = TRUE)
p1=ggplot(data=SV_result,aes(x = Species, y = total_count, group = SV_type, fill = SV_type)) +
  scale_fill_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D"),
    values = c("#e89b27","#00868B","#9DB4CE","#f73e3B","#F9C08A","#A4CB9E"),
    breaks=c("INS","DEL","DUP","CNV","INV","TRANS"),
    labels=c("INS","DEL","DUP","CNV","INV","TRANS"))+
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 100000),    # 设置 X 轴范围
                     breaks = seq(0, 100000, by = 20000), labels = c("0","20k", "40k", "60k", "80k", "100k")) +  # 设置分段
  labs(x = "Species",y = "Number of SV") +coord_flip() +  ##转换横纵坐标
  theme(
    panel.background = element_blank(),
    axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.5, unit = "cm"))
ggsave(p1,filename = "SV_frequency_13genomes.pdf",width = 8,height = 6,dpi=300)
