setwd("F:\\large_data\\gerp_value_B_ZM4\\gerp_sift_paired_info\\gerp_sift_overlap_del_calculate")
SNP=read.table("R_use_format.vcf",header = T,sep = "\t")
group=read.table("D:\\Population_643\\all_761_info\\sift_deleterious_variance\\ind_group.txt",header = T,sep = "\t") ##两列，一列个体名（跟vcf对应），一列分组

##all deleterious variance of individual
mydata=SNP[,3:ncol(SNP)]
het_del=colSums(mydata == "1",na.rm = T)
hom_del=colSums(mydata == "2",na.rm = T)
het_del=as.data.frame(het_del)
hom_del=as.data.frame(hom_del)
combine=cbind(het_del,hom_del)
order=match(rownames(combine),group$taxa)
combine2=cbind(combine,group[order,])
combine2$all_del=combine2$het_del+2*combine2$hom_del
aggregate(combine2$het_del, list(combine2$group), FUN=mean)
aggregate(combine2$hom_del, list(combine2$group), FUN=mean)
aggregate(combine2$all_del, list(combine2$group), FUN=mean)
library(ggpubr)
#t-test
stat.test <- compare_means(
  all_del ~ group, data = combine2,  method = "t.test"
)
stat.test
library(ggplot2)
combine3=combine2[combine2$group=="M.sativa" | combine2$group=="M.caerulea",]
##all info ,outgroup, introgression species
p1<-ggplot(combine3, aes(x=group, y=all_del, fill=group))+labs(x="Species", y="No. of deleterious \n alleles/accession") +
  geom_boxplot(width=0.8)+#ylim(0.7,1.2)+
  scale_fill_manual(
    breaks = c("M.caerulea","M.sativa" ),
    values = c("#9DB4CE", "#F9C08A")
  ) +  ## 不同类型基因对应的颜色
  theme_classic()+theme(
    panel.background = element_blank(),
    legend.position="none",
    axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))
ggsave(plot=p1,"M.sativa_caerulea_deleterious.pdf", width = 6,height = 6,dpi=300)

##bar plot
bar_data=read.table("mean_delterious_variance_per_gene.txt",header = T,sep = "\t")
bar_data$Type2=factor(bar_data$Type2,levels = c("Onefold_essential", "One_copy", "Two_copy", "Three_copy","Four_copy","Fourfold_core"),ordered = TRUE)
bar_data$species=factor(bar_data$species,levels = c("M.caerulea","M.sativa",ordered = TRUE))
p <- ggplot(bar_data, aes(x = Type2, y = mean_del, fill = species)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    breaks = c("M.caerulea","M.sativa" ),
    values = c("#9DB4CE", "#F9C08A")
  ) +  ## 不同类型基因对应的颜色
  labs(title = "", x = "Gene conservation and copy number", y = "Mean number of deleterious alleles/gene") +
  theme_classic()+theme(
    panel.background = element_blank(),
    legend.position="right",
    axis.title=element_text(size=10,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=10,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))
ggsave(plot=p,"M.sativa_caerulea_deleterious_gene_barplot.pdf", width = 8,height = 6,dpi=300)

