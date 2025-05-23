setwd("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\top10_percentage_gene_info/")
library(dplyr)
library(purrr)  # 
file_list <- list.files(path = "./", pattern = "\\.csv", full.names = TRUE)
# 批量读取所有 .txt 文件，将它们存储在一个列表中，并选择所需的列
data_list <- lapply(file_list, function(file) {
  read.csv(file, header = TRUE,sep=",",
           colClasses = c("double",  "character", "character", "character",
                          "character", "character", "character", "character")
           ) 
    
})
merged_top_info <- bind_rows(data_list)
##统计基因出现次数及比例
all_condition_exp <- merged_top_info %>%
  group_by(gene) %>%
  dplyr::summarise(
    count = n(),  # 统计每个类别的出现次数
    Proportion = n()/124 # 计算每个类别的比例,共124种表达条件
  )
write.table(all_condition_exp,"all_condition_exp_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F) ##use to generate manhattan plot

###四个单倍型中不同基因的比例（是否有优势单倍型）
gene_pair=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\combine_four_hap_gene1234_paired_gene.txt",header = T,sep = "\t")
gene_pair_all_top_exp=gene_pair %>% inner_join(all_condition_exp, by=c('gene_name'='gene')) ##提取4个拷贝都有的基因
gene_pair_all_top_exp=gene_pair_all_top_exp[!duplicated(gene_pair_all_top_exp$gene_name),]
result <- gene_pair_all_top_exp %>%
  group_by(hap) %>%
  dplyr::summarise(
    Count = n(),  # 统计每个类别的出现次数
    Proportion = n() / nrow(gene_pair_all_top_exp)  # 计算每个类别的比例
  )
result  ##共17531基因出现在每个条件的前10%
##figures
result$hap=as.factor(result$hap)
p<- ggplot(result, aes(x=hap, y=Proportion, fill=hap)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(
    breaks = c("1", "2", "3", "4"),
    values = c("#F9C08A", "#A4CB9E","#cee2f5","#EDA1A4")
  )+ 
  labs(x = "Haps", y = "Percentage of highly expressed genes") +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype="dashed", color = "blue", size=1.5)+
  theme(
    legend.position = "none",
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p,filename = "proportion_high_exp_gene_top0.1_fpkm.pdf",width = 6,height = 6,dpi=300)

##统计在90%情况下都是高表达的基因
filter_gene=all_condition_exp[all_condition_exp$count>=(length(data_list)*0.9),] ##203个高表达基因
write.table(filter_gene,"203gene_90%condition_high_exp_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)

##基因细节信息
filter_gene_info=merged_top_info[merged_top_info$gene %in% filter_gene$gene,]


##候选基因表达信息查看,Msa017973(优先，激素),Msa036269,Msa040511,Msa048952,Msa091624(光合作用),Msa118329(Zinc finger CCCH domain-containing protein)
##Msa119775(MYB),Msa140769(LRR),Msa154990(cyp450),Msa189697(peroxidase),Msa196066(peroxidase,多种作用)
candidate_gene=filter_gene_info[filter_gene_info$gene=="Msa093079",]
View(candidate_gene)
new_gene=read.table("new_gene.txt",header = T,sep = "\t")
candidate_gene=filter_gene_info[filter_gene_info$gene %in% new_gene$gene,]
result <- candidate_gene %>%
  group_by(gene) %>%
  summarize(mean_value = mean(fpkm))
View(result)


###不同类型基因比例,fourcopy core, onecopy, two copy
group_new=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",head = T,sep = "\t")
filter_gene_info=filter_gene_info %>% inner_join(group_new, by=c('gene'='gene'))

gene_unique <- filter_gene_info[!duplicated(filter_gene_info$gene), ]
result <- gene_unique %>%
  group_by(Type2.y) %>%
  dplyr::summarise(
    Count = n(),  # 统计每个类别的出现次数
    Proportion = n() / nrow(gene_unique)  # 计算每个类别的比例
  )
result
##four copy info
gene_pair=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\combine_four_hap_gene1234_paired_gene.txt",header = T,sep = "\t")
gene_unique2=gene_unique[gene_unique$Type2.x=="Fourfold_core",]
gene_pair_top=gene_pair %>% inner_join(gene_unique2, by=c('gene_name'='gene'))
gene_pair_top=gene_pair_top[!duplicated(gene_pair_top$gene_name),]
write.table(result,"81genes_90%condition_exp_fourfoldgene_detail.txt",row.names = F,col.names = T,sep = "\t",quote = F)
###四个单倍型中高表达fourfold_core基因的比例
result <- gene_pair_top %>%
  group_by(hap) %>%
  dplyr::summarise(
    Count = n(),  # 统计每个类别的出现次数
    Proportion = n() / nrow(gene_pair_top)  # 计算每个类别的比例
  )
result
write.table(result,"90%condition_exp_fourfoldcoregene_four_hap_percent.txt",row.names = F,col.names = T,sep = "\t",quote = F)

##核心基因在不同单倍型表达情况，只在一个单倍型表达还是多个单倍型都表达？
result_exp_hap <- gene_pair_top %>%
  group_by(pair) %>%  ##根据成对信息进行分组
  dplyr::summarise(
    Count = n(),  # 统计每个类别的出现次数
    Proportion = n() / nrow(gene_pair_top)  # 计算每个类别的比例
  )
##81个基因根据pair信息可以归为66个基因（合并pair）
gene_info=gene_pair_top  %>% inner_join(result_exp_hap, by=c('pair'='pair'))
gene_info2=gene_info[gene_info$Count==3,] ##3个hap都表达的基因是哪些
key_gene_exp_info=filter_gene_info[filter_gene_info$gene %in% gene_info2$gene_info,]

result2<- result_exp_hap %>%
  group_by(Count) %>%
  dplyr::summarise(
    number = n(),  # 在几个hap中表达的次数
    Proportion = n() / nrow(result_exp_hap)  # 计算每个类别的比例
  )
result2 #66个基因分别具有1,2,3份拷贝的比例,用于绘制pie fig
new_row <- data.frame(Count = 4,number = 0,Proportion = 0) ##gene expressed in four haps is 0

data1=rbind(result2,new_row)
data1$Count=as.factor(data1$Count)
pie<- ggplot(data1, aes(x="", y=as.numeric(Proportion), fill=Count))+#set basic figure，y信息每个物种都要做
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ #change to pie plot
  scale_fill_manual(values= c("#F9C08A", "#A4CB9E","#cee2f5","#EDA1A4"),name = "Expressed haplotype number",
                    labels = c("1","2","3","4"))+ #change different color for each piece,another color:Zissou1
  labs(x="",y="Expressed haplotype percent")+
  theme_void(base_size = 10, #set size of text
             base_family = "",
             base_line_size = base_size/20,
             base_rect_size = base_size/20)+
  geom_text(aes(x=1.3,size=10,label = paste(format(round(as.numeric(Proportion)*100,2),nsmall=2),"%",sep="")),color = c("black"),position = position_stack(vjust = 0.5),show.legend = FALSE)+ #add text number to pie plot，更换物种记得修改M.group
  theme(plot.margin = margin(t=0, r=0, b=1, l=0, "cm"),
        axis.title = element_text(color="black",size=15),
        axis.text = element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
        legend.title = element_text(color = "black", size = 15), #设置图例标题大小为15，黑色
        legend.text = element_text(color = "black", size = 15),
        legend.position = "right")
ggsave(plot=pie,"Pie_hap1-4_exp_hap.pdf" ,width=10,height=6)

##81个全部条件高表达fourfold core genes在根和叶片中四个单倍型趋势一致？还是不一致？
all_conditions_exp=filter_gene_info[filter_gene_info$gene %in% gene_pair_top$gene_name,]
all_conditions=read.table("124_conditions_info.txt",header = T,sep = "\t")
result<- all_conditions %>%
  group_by(Tissue) %>%
  dplyr::summarise(
    number = n(),  # 统计每个类别的出现次数
    Proportion = n() / nrow(all_conditions)  # 计算每个类别的比例
  )
result ##由于叶片和根的比例占92%，所以后面分析仅用叶片和根组织信息
merge_exp_condition=all_conditions_exp  %>% inner_join(all_conditions, by=c('group'='group',"group2"="group2"))

merge_exp_condition2=merge_exp_condition %>% inner_join(gene_pair_top,by=c('gene'='gene_name'))
result <- merge_exp_condition2 %>%
  group_by(Tissue, hap) %>%
  dplyr::summarise(
    number = n()) %>%  # 统计每个类别的出现次数
  mutate( Proportion = number / sum(number)  # 计算每个 Tissue 内部的比例
  )
result ##叶片，根部的四个hap表达比例，导出到文件，并进行绘图
##figure
##不同基因在不同单倍型中的表达比例
data1= result %>% filter(Tissue %in% c("leaf", "root"))  ##仅使用leaf和root结果，因为别的组织样本量太少
library(ggplot2)
data1$hap=as.factor(data1$hap)
p<- ggplot(data1, aes(x=Tissue, y=Proportion, fill=hap)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(
    breaks = c("1", "2", "3", "4"),
    values = c("#F9C08A", "#A4CB9E","#cee2f5","#EDA1A4")
  )+ 
  labs(x = "Tissue", y = "Percentage of highly expressed genes") +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype="dashed", color = "blue", size=1.5)+
  theme(
    legend.position = "right",
    #panel.background = element_blank(),
    axis.title = element_text(size = 14, colour = 'black'), # 坐标轴标题
    axis.text = element_text(size = 14, colour = 'black'), # 坐标轴标签
    axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )
ggsave(p,filename = "proportion_high_exp_gene_leaf_root.pdf",width = 6,height = 6,dpi=300)


##81个高表达基因在不同条件下（叶片，根部）表达量四个hap间是否具有差异？
result <- merge_exp_condition2 %>%
  group_by(Tissue, hap) %>%
  summarise(count=n(),
    mean_value = mean(fpkm.x, na.rm = TRUE))
result ##不同组织，不同hap的FPKM值差异
##boxplot
merge_exp_condition3=merge_exp_condition2[merge_exp_condition2$Tissue=="leaf" |merge_exp_condition2$Tissue=="root", ]
result <- merge_exp_condition3 %>%
  group_by(gene,Tissue,hap) %>%
  summarise(mean_value = mean(fpkm.x, na.rm = TRUE))
result

library(ggpubr)
stat.test <- compare_means(
  mean_value ~ Tissue, data = result,  method = "t.test"
)
stat.test ##根和叶片直接fpkm是否有差异
stat.test <- compare_means(
  mean_value ~ hap, data = result,  method = "t.test"
)
stat.test ##4个hap之间fpkm是否有显著差异

library(ggplot2)
result$hap=as.factor(result$hap)
p1<-ggplot(result, aes(x=Tissue, y=mean_value, color=hap))+labs(x="Tissue", y="FPKM") +
  geom_boxplot(width=0.8)+ylim(0,1200)+
  scale_color_manual(
    breaks = c("1", "2", "3", "4"),
    values = c("#F9C08A", "#A4CB9E","#cee2f5","#EDA1A4")
  )+ 
  theme_classic()+theme(
    panel.background = element_blank(),
    legend.position="right",
    axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))
ggsave(plot=p1,"leaf_root_exp_four_hap_fpkm.pdf", width = 6,height = 6,dpi=300)
