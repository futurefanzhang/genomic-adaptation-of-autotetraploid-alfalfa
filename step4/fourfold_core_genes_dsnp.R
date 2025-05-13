setwd("F:\\large_data\\gerp_value_B_ZM4\\ZM4_V2_gene_deleterious_info/")

###calculate deleterious genes detail info
gene_pair=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\combine_four_hap_gene1234_paired_gene.txt",header = T,sep = "\t")
##remove redundant mRNA pair (math gene by first transcript)
gene_pair_filtered <-gene_pair %>% filter(mRNA != 1) #extract non first transcript
gene_pair_filtered2 <- gene_pair_filtered %>% group_by(pair) %>% filter(n() == 4) %>% ungroup() ##extract non first transcript in four haps
redundant_gene_pair=unique(gene_pair_filtered2$pair) #redundant gene pair
gene_pair <- gene_pair %>% filter(!(pair %in% redundant_gene_pair)) ##filtered gene_pair
###interact with deleterious SNP for each gene
file_pattern <- "Deleterious_variance_number_Hap.*\\.txt" ##promoter and gene deleterious SNP
files <- list.files(pattern = file_pattern)
df_list <- lapply(files, function(file) {
  read.table(file, header = T, sep = "\t", stringsAsFactors = FALSE)
})
merged_df <- bind_rows(df_list)
mean(merged_df$density)  ##全基因组平均每个基因的有害变异个数是31.64965
quantile(merged_df$density, probs = c(0,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.98,1))
##match tetra_copy core gene
copy_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",header = T,sep = "\t")
all_core_gene_del=merged_df %>% inner_join(copy_info, by=c('gene'='gene'))
result <- all_core_gene_del %>%
  group_by(Type2) %>%
  summarise(total_density = sum(density)) %>%
  mutate(proportion = total_density / sum(total_density))
result
sum(result$total_density)
filtered_df <- merged_df[merged_df$density != 0, ]  ##由于0无法判断是否包括GERP信息，所以不进行计算
mean(filtered_df$density)  ##全基因组非0基因的有害变异平均个数是65.23935

gene_pair_hap1234=gene_pair %>% inner_join(merged_df, by=c('gene_name'='gene')) ##four paired gene deleterious variants
gene_pair_hap1=gene_pair_hap1234[gene_pair_hap1234$hap=="1",]
gene_pair_hap2=gene_pair_hap1234[gene_pair_hap1234$hap=="2",]
gene_pair_hap3=gene_pair_hap1234[gene_pair_hap1234$hap=="3",]
gene_pair_hap4=gene_pair_hap1234[gene_pair_hap1234$hap=="4",]
gene_pair_hap12=gene_pair_hap1 %>% inner_join(gene_pair_hap2, by=c('pair'='pair'))
gene_pair_hap123=gene_pair_hap12 %>% inner_join(gene_pair_hap3, by=c('pair'='pair'))
gene_pair_hap1234_del=gene_pair_hap123 %>% inner_join(gene_pair_hap4, by=c('pair'='pair'))
write.table(gene_pair_hap1234_del,"Deleterious_variance_number_paired_Hap1-4_prometer_gene.txt", row.names = F,quote=F,sep = "\t")
##sum（1-4个基因有害）,汇总成对基因超过平均有害变异位点的个数
rownames(gene_pair_hap1234_del)=gene_pair_hap1234_del$pair
gene_pair_dsnp=gene_pair_hap1234_del[,c("density.x","density.y","density.x.x","density.y.y")]
gene_pair_dsnp$del_hap_count <- apply(gene_pair_dsnp, 1, function(row) {
  # 统计超过65.23935的数值列的数量
  sum(as.numeric(row) > 65.23935, na.rm = TRUE)
})
##match gene pair info and fourfold core info
gene_pair_dsnp$pair=as.integer(rownames(gene_pair_dsnp))
gene_pair_hap1234_four_haps=gene_pair_dsnp %>% inner_join(gene_pair_hap1234_del, by=c('pair'='pair')) ##gene pair info
gene_pair_hap1234_four_haps=gene_pair_hap1234_four_haps %>% inner_join(copy_info, by=c('gene_name.x'='gene')) ##core gene info
fourfold_core_gene_pair_hap1234_four_haps=gene_pair_hap1234_four_haps[gene_pair_hap1234_four_haps$Type2=="Fourfold_core",]
write.table(fourfold_core_gene_pair_hap1234_four_haps,"Fourfold_core_deleterious_variants_paired_Hap1-4_del_prometer_gene.txt", row.names = T,quote=F,sep = "\t")
##export information of fourfold core genes with four hap deleterious variants
gene_pair_hap1234_four_del=fourfold_core_gene_pair_hap1234_four_haps[fourfold_core_gene_pair_hap1234_four_haps$del_hap_count=="4",]
write.table(gene_pair_hap1234_four_del,"Fourfold_core_genes_fourhap_deleterious_variants_number_paired_Hap1-4_four_del_prometer_gene.txt", row.names = T,quote=F,sep = "\t")

##check deleterious haps in fourfold core genes
del_haps_stats <- fourfold_core_gene_pair_hap1234_four_haps %>% group_by(del_hap_count) %>%  # 按 del_haps 分组
  summarise(
    count = n(),  # 计算每个数字的出现次数
    proportion = count / nrow(fourfold_core_gene_pair_hap1234_four_haps)  # 计算出现比例
  )
write.table(del_haps_stats,"Pie_deleterious_hap1-4_percent.txt",sep = "\t",row.names = F)
##figures
del_haps_stats$del_hap_count=as.factor(del_haps_stats$del_hap_count)
pie<- ggplot(del_haps_stats, aes(x="", y=as.numeric(proportion), fill=del_hap_count))+#set basic figure，y信息每个物种都要做
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ #change to pie plot
  scale_fill_manual(values= c("#9DB4CE", "#F9C08A", "#A4CB9E","#cee2f5","#EDA1A4"),name = "Deleterious haplotype number",
                    labels = c("0", "1","2","3","4"))+ #change different color for each piece,another color:Zissou1
  labs(x="",y="Deleterious haplotype percent")+
  theme_void(base_size = 10, #set size of text
             base_family = "",
             base_line_size = base_size/20,
             base_rect_size = base_size/20)+
  geom_text(aes(x=1.3,size=10,label = paste(format(round(as.numeric(proportion)*100,2),nsmall=2),"%",sep="")),color = c("black"),position = position_stack(vjust = 0.5),show.legend = FALSE)+ #add text number to pie plot，更换物种记得修改M.group
  theme(plot.margin = margin(t=0, r=0, b=1, l=0, "cm"),
        axis.title = element_text(color="black",size=15),
        axis.text = element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
        legend.title = element_text(color = "black", size = 15), #设置图例标题大小为15，黑色
        legend.text = element_text(color = "black", size = 15),
        legend.position = "right")
ggsave(plot=pie,"Pie_hap1-4_dele_hap.pdf" ,width=10,height=6)
