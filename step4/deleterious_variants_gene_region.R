setwd("F:\\large_data\\gerp_value_B_ZM4\\ZM4_V2_summary/")
copy_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",header = T,sep = "\t")
detach("package:dplyr", unload = TRUE)
library(dplyr)

##take 3'UTR as example, also need to do for promoter, 5'UTR,CDS,intro
region=read.table("ZM4_V2_intro_region_change_chr.txt",header = F,sep="")
colnames(region)=c("chr", "type", "start", "end", "gene")

chr_list=read.table("chr_list.txt",header = T,sep="\t")
result=list()

for (a in chr_list$chr_type2){
mydata=read.table(paste("F:\\large_data\\gerp_value_B_ZM4\\ZM4_V2_summary\\all_gerp_morethan3\\ZM4_V2",a,"use_gerp_more_than3.txt",sep="_"),header = F,sep="\t") 
colnames(mydata)=c("Position","Neutral_rate","RS_score")
region1=region[region$chr==a,]
print(a)
# 使用两个指针
occ_idx <- 1
num_occurrences <- nrow(mydata)
for (i in 1:nrow(region1)) {
  while (occ_idx <= num_occurrences && mydata$Position[occ_idx] < region1$start[i]) {
    occ_idx <- occ_idx + 1
  }
  occ_start <- occ_idx
  while (occ_idx <= num_occurrences && mydata$Position[occ_idx] <= region1$end[i]) {
    occ_idx <- occ_idx + 1
  }
  region1$gerp_count[i] <- occ_idx - occ_start
}
region1$length=region1$end-region1$start+1
region1$density=(region1$gerp_count/region1$length)*1000 ##the deleterious variance number among every 1k gene
region_info <- region1 %>%
  group_by(gene) %>%
  summarize(
    mean_value = mean(density, na.rm = TRUE))
colnames(region_info)=c("gene","density")

combine_info=region_info %>% inner_join(copy_info, by=c('gene'='gene'))
result1 <- combine_info %>%
  group_by(Type2) %>%
  summarize(
    mean_value = mean(density, na.rm = TRUE),
    sd_value = sd(density, na.rm = TRUE))
result1$chr=a
result[[a]]=result1
}

merged_df <- do.call(rbind, result)
merged_df$region='intro'
write.table(merged_df,"ZM4_V2_intro_deleterious_variance_count.txt",row.names = F,col.names = T,sep = "\t",quote = F)
###merge promoter, 3'UTR, 5'UTR,CDS and intro
promoter=read.table("ZM4_V2_Promoter_deleterious_variance_count.txt",head = T,sep = "\t")
three_prime_UTR=read.table("ZM4_V2_three_prime_UTR_deleterious_variance_count.txt",head = T,sep = "\t")
five_prime_UTR=read.table("ZM4_V2_five_prime_UTR_deleterious_variance_count.txt",head = T,sep = "\t")
CDS=read.table("ZM4_V2_CDS_deleterious_variance_count.txt",head = T,sep = "\t")
intro=read.table("ZM4_V2_intro_deleterious_variance_count.txt",head = T,sep = "\t")
merged_region=rbind(promoter,three_prime_UTR,five_prime_UTR,CDS,intro)
write.table(merged_region,"summary_deleterious_different_region.txt",row.names = F,col.names = T,sep = "\t",quote = F)

###summary figure
pan_core=read.table("summary_deleterious_different_region.txt",head = T,sep = "\t")
result1 <- pan_core %>%
  group_by(Type2,region) %>%
  summarize(
    mean_value = mean(mean_value, na.rm = TRUE),
    proportion = mean_value/ sum(mean_value)
    ) ##generate the mean deleterious SNP number of each chr
result1 
result <- result1 %>%
  group_by(Type2) %>%
  mutate(total_mean = sum(mean_value)) %>%
  ungroup() %>%
  mutate(proportion = mean_value / total_mean) ##calculate the proportion of deleterious SNP in different regions
write.table(result,"summary_each_type_gene_deleterious_variance_count_proportion.txt",row.names = F,col.names = T,sep = "\t",quote = F)
result <- result1 %>%
  group_by(Type2) %>%
  summarize(total_mean = sum(mean_value)) %>%
  mutate(proportion = total_mean / sum(total_mean))  ##chick which group (onefold essential or fourfold core) genes have the most deleterious SNP

library(ggplot2)
result1$Type2 <- str_replace_all(result1$Type2, c("one_copy" = "One_copy", "two_copy" = "Two_copy","three_copy" = "Three_copy","four_copy" = "Four_copy"))
result1$Type2 <- factor(result1$Type2,levels = c("Onefold_essential","One_copy","Two_copy","Three_copy","Four_copy","Fourfold_core"),ordered = TRUE) 
result1$region <- str_replace_all(result1$region, c("intro" = "Intron"))
result1$region <- factor(result1$region,levels = c("Promoter","5-Primer UTR","CDS","Intron","3-Primer UTR"),ordered = TRUE)
p1=ggplot(data=result1,aes(x = Type2, y = mean_value, group = region, fill = region)) +
  scale_fill_manual(
    breaks = c("Promoter","5-Primer UTR","CDS","Intron","3-Primer UTR"),
    values = c("#9DB4CE", "#F9C08A", "#EDA1A4", "#A4CB9E","#cee2f5")
  ) +  geom_bar(stat = "identity") +
  labs(x = "",y = "Number of deleterious variants per Kb") +
  theme(
    panel.background = element_blank(),
    axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.5, unit = "cm"))
ggsave(p1,filename = "Number of deleterious variance_core_gene_region2.pdf",width = 12,height = 6,dpi=300)  ##Figure 4B
