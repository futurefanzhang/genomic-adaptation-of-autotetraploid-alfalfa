#use the results of orthofinders
##1: /data/home/zhangfan/compare_genomoe_medicago/genespace_new/orthofinder/Results_May09/Orthogroups/Orthogroups.GeneCount.tsv
##2: /data/home/zhangfan/compare_genomoe_medicago/genespace_new/orthofinder/Results_May09/Comparative_Genomics_Statistics/Statistics_Overall.tsv
##analysis in R, four figures, barplot, pieplot, heatmap, and Growth curve figures
library(ggplot2)
setwd("E:\\Genome\\Comparative_genomics\\core_dispensable_genes_v2/")
dat = read.csv("bar_plot.txt",header = T,sep = "\t")
###图1：柱形图
#bar plot
library(ggplot2)
dat$Species <- factor(dat$Species,levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13"),ordered = TRUE)
dat$Group <- factor(dat$Group,levels = c("Private", "Dispensable","Softcore", "Core"),ordered = TRUE)
p <- ggplot(dat,aes(y=gene_family,x=Species,fill=Group)) +
  geom_bar(stat="identity",position = "dodge",width=0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,15000),breaks = c(0,5000,10000,15000))+
                                                                      
  labs(x = "Genome number",y = "Number of gene families") +
  scale_fill_manual(breaks=c("Private", "Dispensable","Softcore", "Core"),
    values=c("#9DB4CE","#F9C08A","#EDA1A4","#A4CB9E"))+
  theme(
        panel.background = element_blank(),
        axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
        axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(size=0.5, colour = 'black'), #轴线
        plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.5, unit = "cm")) 

ggsave(p,filename = "bar.pdf",width = 9,height = 6,dpi=300)
#ggsave(p,filename = "bar.jpg",width = 9,height =6,dpi=300)
###图2：饼图
#pie plot
data1 = read.csv("pie_plot.txt",header = T,sep = "\t")
data1$Group <- factor(data1$Group,levels = c("Private", "Dispensable","Softcore", "Core"),ordered = TRUE)
pie<- ggplot(data1, aes(x="", y=Percentage, fill=Group))+#set basic figure
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ #change to pie plot
  scale_fill_manual(breaks=c("Private", "Dispensable","Softcore", "Core"),
                    values=c("#9DB4CE","#F9C08A","#EDA1A4","#A4CB9E"))+
   theme_void(base_size = 20, #set size of text
             base_family = "",
             base_line_size = base_size/20,
             base_rect_size = base_size/20)+
  geom_text(aes(x=1.2,label = paste(format(round(Percentage,2),nsmall=2),"%",sep="")),size=8,color = c("black"),position = position_stack(vjust = 0.5))+ #add text number to pie plot
  theme(plot.margin = unit(x = c(0, 4, 0, 0), units = "mm"),#change the position of legend
        legend.title = element_text(color = "black", size = 20),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(color = "black", size = 20),
        legend.key.width=unit(0.5,"cm"))#change legend size
       
ggsave(pie,filename = "pie_plot.pdf",width = 9,height = 6,dpi=300)

###图3：热图
##distribution heatmap plot
mydata=read.table("Orthogroups.GeneCount.tsv",header = T,sep = "\t")
rownames(mydata)=mydata$Orthogroup
mydata <- subset(mydata, select = -c(Orthogroup,Total))  ##remove: Orthogroup,Total
mydata[] <- sapply(mydata, function(x) ifelse(x > 1, 1, x)) ##修改超过1的数值为1
mydata$sum=rowSums(mydata)
mydata <- mydata[order(-mydata$sum),]
mydata$sum[mydata$sum < 11 & mydata$sum > 1] <- "Dispensable"
mydata$sum[mydata$sum == 13] <- "Core"
mydata$sum[mydata$sum == 1] <- "Private"
mydata$sum[mydata$sum == 11 | mydata$sum == 12] <- "Softcore"
write.table(mydata,"pan_gene_info.txt",row.names = T,col.names = T,sep = "\t",quote = F) ##导出处理后的文件
rowname_use=as.data.frame(mydata$sum)
rowname_use$number=seq(1:nrow(rowname_use))
colnames(rowname_use)=c("sum","number")
library(dplyr)
x_labinfo <- rowname_use %>%
  group_by(sum) %>%
  mutate(count = n()) %>%
  slice(1) %>%
  ungroup()  # 用于移除分组，使后续的操作不受分组影响
x_labinfo$lab=x_labinfo$count+x_labinfo$number-1
colnames(x_labinfo)=c("Gene","start","count","end")
x_lab1=x_labinfo[,c("Gene","start")]
x_lab2=x_labinfo[,c("Gene","end")]
colnames(x_lab2)=c("Gene","start")
x_lab=rbind(x_lab1,x_lab2)
mydata=subset(mydata, select = -c(sum))  
mydata$RowName <- seq_len(nrow(mydata))

library(reshape2)
mydata_long <- melt(mydata, id.vars = "RowName")
colnames(mydata_long) <- c("RowName", "ColumnName", "Value")
color_range <- colorRampPalette(c("#FF7E4F","#1C9FE2"))(2)
#color_range <- colorRampPalette(c("#F8766D","#7CAE00"))(2)
mydata_long$Value=factor(mydata_long$Value,levels = c("1", "0"),ordered = TRUE)
library(ggplot2)

p1<-ggplot(mydata_long, aes(x = RowName, y = ColumnName,fill=Value)) +
  geom_tile() +
  scale_fill_manual(values = color_range,name = "",labels = c("Present", "Absent")) +
  scale_x_continuous(expand = expansion(mult = 0), 
    breaks = x_lab$start,  # 位置
    labels =  x_lab$Gene  # 标签
  )+
  labs(x = "", y = "") +
  theme(
    legend.position = "top",
    panel.background = element_blank(),
    axis.text.x=element_blank(),
    #axis.text.x = element_text(size = 15, color = "Black",hjust = -1,vjust=4),
    axis.text.y = element_text(size = 15,color = "Black"),
    axis.title.x=element_blank(),
    axis.ticks.x = element_line(size = 0.5),  # 调整 X 轴刻度线的粗细
    axis.ticks.y = element_line(size = 0.5),  # 调整 Y 轴刻度线的粗细
    axis.ticks.length = unit(0.4, "cm") 
    #axis.ticks.x=element_blank()
  )
# 保存图像,需要AI调整一下重叠字，不同基因组排列顺序：提前把导入的表格按照自己希望的顺序排列
ggsave(plot=p1,"pan_genome_heatmap.pdf", width = 10, height = 6, dpi = 300)

###图4：饱和曲线图
###pan and core gene point figure
mydata=read.table("pan_gene_info.txt",row.names = 1,head = T,sep = "\t")
df=subset(mydata, select = -c(sum))
n_cols <- ncol(df)

# 创建一个函数，用于从数据框中提取指定列数的列组合，并创建新的数据框
create_dfs <- function(cols) {
  # 使用combn生成所有可能的列组合
  column_combinations <- combn(names(df), cols)
  
  # 使用lapply遍历每个组合，并为每个组合创建一个新的数据框
  new_dfs <- lapply(seq_along(column_combinations[1, ]), function(i) {
    data.frame(
      df[column_combinations[, i]]
    )
  })
  # 返回新创建的数据框列表
  
  return(new_dfs)
}

# 创建一个列表，用于存储所有新创建的数据框
all_new_dfs <- list()

# 遍历从1到数据框列数的所有列数，创建相应的新数据框，这一步费时间，我的是13个基因组，如果更多基因组，会耗时久
for (cols in 1:n_cols) {
  all_new_dfs[[as.character(cols)]] <- create_dfs(cols)
}

##计算出core和pan基因情况,把上面的列表转换为数据框，然后统计每种组合的core和pan基因情况
use_list=list()
for (a in 1:length(all_new_dfs)){
length1=all_new_dfs[[a]]
use_dataframe=data.frame()
for (i in 1:length(length1)){
dataframe=as.data.frame(length1[[i]])
dataframe$sum=rowSums(dataframe)
core=sum(dataframe$sum==(ncol(dataframe)-1))  ##所有个体都出现，是core基因
pan=sum(dataframe$sum!=0)   ##有一个个体出现，pan基因
use_dataframe[i,1]=i
use_dataframe[i,2]=core
use_dataframe[i,3]=pan
}
use_list[[a]]=use_dataframe
}
# 提取子集的编号
subset_numbers <- as.character(1:length(use_list))

# 创建一个函数，用于将子集转换为数据框，并在每个子集中添加一个编号列
convert_to_df <- function(subset, number) {
  subset$SubsetNumber <- number
  return(subset)
}

# 将列表中的子集合并成一个数据框，并为每个子集内部的数据框添加编号列
df <- do.call(rbind, lapply(seq_along(use_list), function(i) {
  convert_to_df(use_list[[i]], subset_numbers[i])
}))
colnames(df)=c("combination_type","core_gene","pan_gene","genome_number")
write.table(df,"pan_and_core_gene_info.txt",row.names = F,col.names = T,sep = "\t",quote = F)
##figure
pan_gen=subset(df, select = c(pan_gene,genome_number))
pan_gen$Type="Pan"
colnames(pan_gen)=c("Gene_number","Genome_number","Type")
core_gen=subset(df, select = c(core_gene,genome_number))
core_gen$Type="Core"
colnames(core_gen)=c("Gene_number","Genome_number","Type")
all_info=rbind(core_gen,pan_gen)
all_info$Genome_number <- factor(all_info$Genome_number,levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13"),ordered = TRUE)

##
p1<-ggplot(data=all_info,aes(x=Genome_number,y=Gene_number,colour=Type))+
  geom_point(size=1)+
  labs(x="Genome number",y="Number of gene families")+#PC1_7.97%,PC2_7.20%,PC3_6.04%
  scale_color_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D","gray","gray"),
    values = c("#FF7E4F","#1C9FE2"),
    breaks=c("Core", "Pan"),
    labels=c("Core", "Pan"))+
  stat_summary(
    aes(group = Type),  # 根据 Type 分组
    fun = mean,          # 计算均值
    geom = "line",       # 绘制线
    linewidth = 1             # 线的粗细
  )+
  #labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=10),
                        axis.title.y = element_text(color="black",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 10), #设置图例标题大小为20，黑色
                        legend.text = element_text(color = "black", size = 10),
                        legend.position = "right")
ggsave(plot=p1,"Core_pan_gene_number_compare.pdf", width = 10, height = 6)
