setwd("E:\\Genome\\PGGB\\perennal_annual_GWAS")
col_types <- c("character", rep("numeric", 19))

chr="chr1"
geno=read.table(paste(chr,"sv_genotype.txt",sep="_"),header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))

ana_data=geno[,3:ncol(geno)]
ana_data=ana_data[,c("Mcalanda","Mcalong","Mt5","MtHM078","R108hic","ZM1","ZM4V22","ZM4V23","ZM4V24","xjdy1","xjdy2","xjdy3","xjdy4")]
ana_data=na.omit(ana_data)
ana_data=t(ana_data)
pca_result <- prcomp(ana_data, center = TRUE, scale. = FALSE)
##PCA得分
# 1. 提取标准差
std_dev <- pca_result$sdev  # 主成分的标准差

# 2. 计算每个主成分解释的方差
variance <- std_dev^2  # 方差

# 3. 计算解释的百分比
explained_variance_ratio <- variance / sum(variance) * 100  # 每个主成分解释的百分比

# 4. 创建数据框以便于查看
explained_variance_df <- data.frame(PC = paste0("PC", 1:length(explained_variance_ratio)),
                                    Variance = variance,
                                    ExplainedVariancePercentage = explained_variance_ratio)
pc1_percentage <- explained_variance_df$ExplainedVariancePercentage[1]
pc2_percentage <- explained_variance_df$ExplainedVariancePercentage[2]
# 6. 可视化 PCA 结果
# 将 PCA 得分转换为数据框
pca_scores <- as.data.frame(pca_result$x)

# 为了便于可视化，将物种信息加入数据框
pca_scores$Species <- rownames(pca_scores)

# 绘制 PCA 图
library(ggplot2)
p1=ggplot(data = pca_scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  labs(title = "",
       x = paste("PC1 (", round(pc1_percentage, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pc2_percentage, 2), "%)", sep = "")) +
  scale_color_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D"),
    values = c("#FABB2E","#19B700","#3ea32e","#00868B","#42d4f4","#BF3EFF","#742b99","#A3A500","#dbde43","#6c6e09","#ee0000","#f73e3B","#fa6b6b","#fc8181","#cf4a4a","#940303","#700202","#750e0e","#661a1a"),
    breaks=c("Mara", "Mcalanda","Mcalong","Mlup1", "Mpoly1","Mrulanda","Mruzhiwusuo","Mt5","MtHM078","R108hic","xjdy1","xjdy2","xjdy3","xjdy4","ZM1","ZM4V22","ZM4V23","ZM4V24"),
    labels=c("Mara", "Mcalanda","Mcalong","Mlup1", "Mpoly1","Mrulanda","Mruzhiwusuo","Mt5","MtHM078","R108hic","xjdy1","xjdy2","xjdy3","xjdy4","ZM1","ZM4V22","ZM4V23","ZM4V24"))+
  labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "right")
ggsave(plot=p1,paste(chr,"SV_PC_1_2_mt_msa.pdf",sep="_"), width=7.5,height=6)

##合并全部染色体，做汇总的PCA
#由于Chr5,chr6和chr7包含的SV变异太大，无法包含全部个体，因此仅用其余染色体
col_types <- c("character", rep("numeric", 19))

geno1=read.table("chr1_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))
geno2=read.table("chr2_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))
geno3=read.table("chr3_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))
geno4=read.table("chr4_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))
geno8=read.table("chr8_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))
col_types <- c("character", rep("numeric", 21))
geno5=read.table("chr5_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))
col_types <- c("character", rep("numeric", 17))
geno6=read.table("chr6_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))
col_types <- c("character", rep("numeric", 18))
geno7=read.table("chr7_sv_genotype.txt",header = T,sep = "\t",colClasses = col_types,na.strings = c(".", "", "NA"))

geno=rbind(geno1,geno2,geno3,geno4,geno8)

ana_data=geno[,3:ncol(geno)] ##all sample
str(ana_data)
ana_data$xjdy4=as.numeric(ana_data$xjdy4)

ana_data[ana_data > 2] <- NA
ana_data=na.omit(ana_data)



ana_data=t(ana_data)

pca_result <- prcomp(ana_data, center = TRUE, scale. = FALSE)
std_dev <- pca_result$sdev  # 主成分的标准差
variance <- std_dev^2  # 方差
explained_variance_ratio <- variance / sum(variance) * 100  # 每个主成分解释的百分比

pc1_percentage <- explained_variance_ratio[1]
pc2_percentage <- explained_variance_ratio[2]
# 6. 可视化 PCA 结果
# 将 PCA 得分转换为数据框
pca_scores <- as.data.frame(pca_result$x)

# 为了便于可视化，将物种信息加入数据框
pca_scores$Species <- rownames(pca_scores)

# 绘制 PCA 图
library(dplyr)
pca_scores <- pca_scores %>% mutate(Species = recode(Species,
 "Mara" = "Mara","Mcalanda" = "Mca_landa", "Mcalong" = "Mca_long","Mlup1" = "Mlup","Mpoly1" = "Mpoly","Mrulanda" = "Mru_landa",
  "Mruzhiwusuo" = "Mru_zhiwusuo","Mt5" = "Mt_A17","MtHM078" = "Mt_HM078","R108hic" = "Mt_R108", "xjdy1" = "XJDY_hap1","xjdy2" = "XJDY_hap2", "xjdy3" = "XJDY_hap3",
"xjdy4" = "XJDY_hap4", "ZM1" = "Zhongmu1", "ZM4V22" = "ZM4_V2_hap2","ZM4V23" = "ZM4_V2_hap3","ZM4V24" = "ZM4_V2_hap4"))
library(ggplot2)
library(ggrepel)
p1=ggplot(data = pca_scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  #geom_text_repel(aes(label = Species),box.padding = unit(0.6, "lines"),point.padding = unit(0.1, "lines"),
  #                size = 3) +
  labs(title = "",
       x = paste("PC1 (", round(pc1_percentage, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pc2_percentage, 2), "%)", sep = "")) +
  scale_color_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D"),
    values = c("#FABB2E","#19B700","#3ea32e","#00868B","#42d4f4","#BF3EFF","#742b99","#A3A500","#dbde43","#6c6e09","#ee0000","#f73e3B","#fa6b6b","#fc8181","#e89b27","#700202","#750e0e","#661a1a"),
    breaks=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2_hap2","ZM4_V2_hap3","ZM4_V2_hap4"),
    labels=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2_hap2","ZM4_V2_hap3","ZM4_V2_hap4"))+
  labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "left")
ggsave(plot=p1,"chr_all_SV_PC_1_2.pdf", width=8,height=6)

##M. sativa-falcata complex
geno=rbind(geno1,geno2,geno3,geno4,geno8)

ana_data=geno[,3:ncol(geno)] ##all sample
ana_data=ana_data[,c("Mcalanda","Mcalong","ZM1","ZM4V22","ZM4V23","ZM4V24","xjdy1","xjdy2","xjdy3","xjdy4")]
str(ana_data)
ana_data$xjdy4=as.numeric(ana_data$xjdy4)
ana_data[ana_data > 2] <- NA
ana_data=na.omit(ana_data)



ana_data=t(ana_data)

pca_result <- prcomp(ana_data, center = TRUE, scale. = FALSE)
std_dev <- pca_result$sdev  # 主成分的标准差
variance <- std_dev^2  # 方差
explained_variance_ratio <- variance / sum(variance) * 100  # 每个主成分解释的百分比

pc1_percentage <- explained_variance_ratio[1]
pc2_percentage <- explained_variance_ratio[2]
# 6. 可视化 PCA 结果
# 将 PCA 得分转换为数据框
pca_scores <- as.data.frame(pca_result$x)

# 为了便于可视化，将物种信息加入数据框
pca_scores$Species <- rownames(pca_scores)

# 绘制 PCA 图
library(dplyr)
pca_scores <- pca_scores %>% mutate(Species = recode(Species,
                                                     "Mara" = "Mara","Mcalanda" = "Mca_landa", "Mcalong" = "Mca_long","Mlup1" = "Mlup","Mpoly1" = "Mpoly","Mrulanda" = "Mru_landa",
                                                     "Mruzhiwusuo" = "Mru_zhiwusuo","Mt5" = "Mt_A17","MtHM078" = "Mt_HM078","R108hic" = "Mt_R108", "xjdy1" = "XJDY_hap1","xjdy2" = "XJDY_hap2", "xjdy3" = "XJDY_hap3",
                                                     "xjdy4" = "XJDY_hap4", "ZM1" = "Zhongmu1", "ZM4V22" = "ZM4_V2_hap2","ZM4V23" = "ZM4_V2_hap3","ZM4V24" = "ZM4_V2_hap4"))
library(ggplot2)
library(ggrepel)
p1=ggplot(data = pca_scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Species),box.padding = unit(0.6, "lines"),point.padding = unit(0.1, "lines"),
                  size = 3) +
  labs(title = "",
       x = paste("PC1 (", round(pc1_percentage, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pc2_percentage, 2), "%)", sep = "")) +
  scale_color_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D"),
    values = c("#FABB2E","#19B700","#3ea32e","#00868B","#42d4f4","#BF3EFF","#742b99","#A3A500","#dbde43","#6c6e09","#ee0000","#f73e3B","#fa6b6b","#fc8181","#e89b27","#700202","#750e0e","#661a1a"),
    breaks=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2_hap2","ZM4_V2_hap3","ZM4_V2_hap4"),
    labels=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2_hap2","ZM4_V2_hap3","ZM4_V2_hap4"))+
  labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "none")
ggsave(plot=p1,"chr_all_SV_PC_1_2_medicago_sativa_falcata_complex.pdf", width=8,height=6)

##single chr
geno=geno1 #chr1
##geno=geno2,chr2
##the same for the rest
ana_data=geno[,3:ncol(geno)] ##all sample
ana_data=ana_data[,c("Mcalanda","Mcalong","ZM1","ZM4V22","ZM4V23","ZM4V24","xjdy1","xjdy2","xjdy3","xjdy4")]
str(ana_data)
ana_data$xjdy4=as.numeric(ana_data$xjdy4)
ana_data[ana_data > 2] <- NA
ana_data=na.omit(ana_data)



ana_data=t(ana_data)

pca_result <- prcomp(ana_data, center = TRUE, scale. = FALSE)
std_dev <- pca_result$sdev  # 主成分的标准差
variance <- std_dev^2  # 方差
explained_variance_ratio <- variance / sum(variance) * 100  # 每个主成分解释的百分比

pc1_percentage <- explained_variance_ratio[1]
pc2_percentage <- explained_variance_ratio[2]
# 6. 可视化 PCA 结果
# 将 PCA 得分转换为数据框
pca_scores <- as.data.frame(pca_result$x)

# 为了便于可视化，将物种信息加入数据框
pca_scores$Species <- rownames(pca_scores)

# 绘制 PCA 图
library(dplyr)
pca_scores <- pca_scores %>% mutate(Species = recode(Species,
                                                     "Mara" = "Mara","Mcalanda" = "Mca_landa", "Mcalong" = "Mca_long","Mlup1" = "Mlup","Mpoly1" = "Mpoly","Mrulanda" = "Mru_landa",
                                                     "Mruzhiwusuo" = "Mru_zhiwusuo","Mt5" = "Mt_A17","MtHM078" = "Mt_HM078","R108hic" = "Mt_R108", "xjdy1" = "XJDY_hap1","xjdy2" = "XJDY_hap2", "xjdy3" = "XJDY_hap3",
                                                     "xjdy4" = "XJDY_hap4", "ZM1" = "Zhongmu1", "ZM4V22" = "ZM4_V2_hap2","ZM4V23" = "ZM4_V2_hap3","ZM4V24" = "ZM4_V2_hap4"))
library(ggplot2)
library(ggrepel)
p1=ggplot(data = pca_scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Species),box.padding = unit(0.6, "lines"),point.padding = unit(0.1, "lines"),
                  size = 3) +
  labs(title = "",
       x = paste("PC1 (", round(pc1_percentage, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pc2_percentage, 2), "%)", sep = "")) +
  scale_color_manual(#values = c("#9590FF","#BFEFFF","#EEEE00","#F8766D","#FABB2E","#A3A500","#00868B","#A8422D"),
    values = c("#FABB2E","#19B700","#3ea32e","#00868B","#42d4f4","#BF3EFF","#742b99","#A3A500","#dbde43","#6c6e09","#ee0000","#f73e3B","#fa6b6b","#fc8181","#e89b27","#700202","#750e0e","#661a1a"),
    breaks=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2_hap2","ZM4_V2_hap3","ZM4_V2_hap4"),
    labels=c("Mara", "Mca_landa","Mca_long","Mlup", "Mpoly","Mru_landa","Mru_zhiwusuo","Mt_A17","Mt_HM078","Mt_R108","XJDY_hap1","XJDY_hap2","XJDY_hap3","XJDY_hap4","Zhongmu1","ZM4_V2_hap2","ZM4_V2_hap3","ZM4_V2_hap4"))+
  labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "none")
ggsave(plot=p1,"chr_all_SV_PC_1_2_medicago_sativa_falcata_complex_chr1.pdf", width=8,height=6)
