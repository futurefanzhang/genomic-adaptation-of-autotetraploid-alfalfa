setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\all_top_percent_gene_proportion")
library(dplyr)
library(purrr)  # 
file_list <- list.files(path = "./", pattern = "\\gene_exp.txt", full.names = TRUE)
treat_group=read.table("124_conditions_info.txt",header = T,sep = "\t")
# 批量读取所有 .txt 文件，将它们存储在一个列表中，并选择所需的列
data_list=list()
data_list <- lapply(file_list, function(file) {
  print(file)  # 打印当前处理的文件名
  # 读取数据并合并
  combined_info_total <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    inner_join(treat_group, by = c('group' = 'group', 'group2' = 'group2'))  # 合并数据
  # 选择特定列
  combined_info_total %>%
    select(Type2, proportion, top_percent, Treatment_group2)  
})
merged_top_info <- bind_rows(data_list)
##figure
merged_top_info$top_percent2=1-merged_top_info$top_percent
merged_top_info$top_percent2=as.factor(merged_top_info$top_percent2)
merged_top_info$top_percent2=factor(merged_top_info$top_percent2,levels = c("0.5", "0.25", "0.1", "0.05","0.01"),ordered = TRUE)
detach("package:dplyr", unload = TRUE)
library(dplyr)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
write.table(merged_top_info,"proportion_accession_toppercent_fpkm.txt", row.names = F,col.names = T,sep = "\t",quote = F) 
df3 <- data_summary(merged_top_info, varname="proportion", groupnames=c("Type2","top_percent2", "Treatment_group2"))
colnames(df3)=c("Type2","top_percent2", "Treatment","proportion","sd")
background=read.table("raw_gene_type2_percent.txt",header = T,sep = "\t")                  
df3=rbind(df3,background)
df3$top_percent2=as.factor(df3$top_percent2)
df3$top_percent2=factor(df3$top_percent2,levels = c("0.5", "0.25", "0.1", "0.05","0.01"),ordered = TRUE)
df3$Type2 <- str_replace_all(df3$Type2, c("one_copy" = "One_copy", "two_copy" = "Two_copy","three_copy" = "Three_copy","four_copy" = "Four_copy"))

df3$Type2=factor(df3$Type2,levels = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),ordered = TRUE)
df3$Treatment=factor(df3$Treatment,levels = c("Biotic stress","Abiotic stress", "Normal","Genomic background"),ordered = TRUE)
write.table(df3,"all_top_exp_percent.txt",row.names = F,quote=F,sep = "\t")
p1=ggplot(df3, aes(x=top_percent2, y=proportion, group=Treatment,color=Treatment)) + 
  scale_color_manual(
    breaks = c("Biotic stress","Abiotic stress", "Normal","Genomic background"),
    values = c("#EE3377","#33BBEE","#009988" ,"gray60")
  ) +  
  geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+facet_grid(. ~ Type2)+
  labs(x = "Top expressed gene percentage", y = "Percentage") +
  #scale_color_brewer(palette="Paired")+
  #theme_classic()+
  theme(
    panel.background = element_rect(fill = 'white', colour = 'lightgray'),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, colour = 'black', angle = 30,vjust = 0.5),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p1,filename = "proportion_accession_toppercent_fpkm3.pdf",width = 15,height = 5,dpi=300)
