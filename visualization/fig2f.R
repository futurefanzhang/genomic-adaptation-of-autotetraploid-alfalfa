setwd("E:\\Genome\\PGGB\\perennal_annual_GWAS")
library(dplyr)
combined_data <- data.frame()
for (i in 1:8){
chr1=read.table(paste0("chr",i,"_sv_genotype.txt"),header = T,sep = "\t",na.strings = c(".", "", "NA"))
xjdy=chr1[,c("POS","xjdy1","xjdy2","xjdy3","xjdy4")]
xjdy$sv_number=4-rowSums(is.na(xjdy))
xjdy$sv_number[xjdy$sv_number==0]<-NA
ZM4_V2=chr1[,c("POS","ZM4V22","ZM4V23","ZM4V24")]
ZM4_V2$sv_number=3-rowSums(is.na(ZM4_V2))
ZM4_V2$sv_number[ZM4_V2$sv_number==0]<-NA
chr1_13genome <- chr1 %>% select(-c(xjdy1, xjdy2, xjdy3, xjdy4, ZM4V22, ZM4V23, ZM4V24))
chr1_13genome=merge(chr1_13genome, xjdy[, c("POS", "sv_number")], by = "POS", all.x = TRUE)
chr1_13genome=chr1_13genome %>% rename(xjdy = sv_number)
chr1_13genome=merge(chr1_13genome, ZM4_V2[, c("POS", "sv_number")], by = "POS", all.x = TRUE)
chr1_13genome=chr1_13genome %>% rename(ZM4_V2 = sv_number)
chr1_13genome_sv_miss=rowSums(is.na(chr1_13genome))
chr1_13genome$present=rowSums(!is.na(chr1_13genome[, 3:ncol(chr1_13genome)])) ##all 13 genomes sv distribution
use_info=chr1_13genome[,c("CHROM","POS","present")]
use_info=use_info[use_info$present>0,]
combined_data=rbind(combined_data,use_info)
}
combined_data_use=as.data.frame(combined_data)
library(dplyr)
frequency_df_all <- combined_data_use %>%
  count(present)
write.table(frequency_df_all,"all_13genome_sv_info.txt",row.names = F,col.names = T,sep = "\t",quote = F)
##medicago sativa complex
combined_data2 <- data.frame()
for (i in 1:8){
chr1=read.table(paste0("chr",i,"_sv_genotype.txt"),header = T,sep = "\t",na.strings = c(".", "", "NA"))
xjdy=chr1[,c("POS","xjdy1","xjdy2","xjdy3","xjdy4")]
xjdy$sv_number=4-rowSums(is.na(xjdy))
xjdy$sv_number[xjdy$sv_number==0]<-NA
ZM4_V2=chr1[,c("POS","ZM4V22","ZM4V23","ZM4V24")]
ZM4_V2$sv_number=3-rowSums(is.na(ZM4_V2))
ZM4_V2$sv_number[ZM4_V2$sv_number==0]<-NA
chr1_13genome <- chr1 %>% select(-c(xjdy1, xjdy2, xjdy3, xjdy4, ZM4V22, ZM4V23, ZM4V24))
chr1_13genome=merge(chr1_13genome, xjdy[, c("POS", "sv_number")], by = "POS", all.x = TRUE)
chr1_13genome=chr1_13genome %>% rename(xjdy = sv_number)
chr1_13genome=merge(chr1_13genome, ZM4_V2[, c("POS", "sv_number")], by = "POS", all.x = TRUE)
chr1_13genome=chr1_13genome %>% rename(ZM4_V2 = sv_number)
chr1_13genome_sv_miss=rowSums(is.na(chr1_13genome))
chr1_13genome$present=rowSums(!is.na(chr1_13genome[, 3:ncol(chr1_13genome)])) ##all 13 genomes sv present
chr1_13genome=chr1_13genome[chr1_13genome$present>0,]
chr1_5genome=chr1_13genome[,c("Mcalanda","Mcalong","ZM1","xjdy","ZM4_V2","present")]
chr1_5genome_sv_miss=rowSums(is.na(chr1_5genome))
chr1_5genome$present2=5-chr1_5genome_sv_miss ## five genomes sv present
filtered_data <- chr1_5genome %>% filter(present == present2) ##extract the present info only present in M.sativa complex
combined_data2=rbind(combined_data2,filtered_data)
}
combined_data_use=as.data.frame(combined_data2)
frequency_df_all2 <- combined_data_use %>%
  count(present2)
write.table(frequency_df_all2,"all_5genome_sv_info.txt",row.names = F,col.names = T,sep = "\t",quote = F)

##combine previous two tables, generate figure
combined_data <- full_join(frequency_df_all, frequency_df_all2, by = c("present" = "present2"))
colnames(combined_data)=c("Frequency","All 13 genomes","M. sativa-falcata complex")
combined_data$`M. sativa-falcata complex`[is.na(combined_data$`M. sativa-falcata complex`)] <- 0
part1=combined_data[,c(1,2)]
colnames(part1)=c("Frequency","Count")
part1$Type="All 13 genomes"
part2=combined_data[,c(1,3)]
colnames(part2)=c("Frequency","Count")
part2$Type="M. sativa-falcata complex"
all_info=rbind(part1,part2)
##figures
all_info$Frequency=as.factor(all_info$Frequency)
write.table(all_info,"all_SV_frequency.txt",row.names = F,col.names = T,sep = "\t",quote = F)

p1=ggplot(data=all_info, aes(x=Frequency, y=Count, fill=Type, color=Type, alpha=Type)) +
  geom_bar(stat="identity", position ="identity") +
  scale_colour_manual(values=c("#027fbf", "#ff4805")) +
  scale_fill_manual(values=c("#1C9FE2", "#ff7e4f")) +
  scale_alpha_manual(values=c(.3, .7))+ 
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "top")
ggsave(plot=p1,"all_SV_frequency.pdf", width=8,height=6)
