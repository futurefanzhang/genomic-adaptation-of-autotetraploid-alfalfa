###core gene info
pan_core_info=read.table("E:\\Genome\\Comparative_genomics\\core_dispensable_genes_v2\\all_gene_orthogroup_pan_gene_info.txt",head = T,sep = "\t")
pan_core_info=pan_core_info %>% separate(ZM4_V2, into = c("gene_name", "mRNA"), sep = ".t")
pan_core_info$sum=factor(pan_core_info$sum,levels = c("Core","Softcore","Dispensable","Private"),ordered = TRUE)

##核对基因保留是否会出错
#unique_variants <- pan_core_info %>% group_by(gene_name) %>% filter(n_distinct(sum) > 1) %>% ungroup() ##查找重复基因
#distinct_rows <- unique_variants %>% group_by(gene_name) %>% distinct(sum, .keep_all = TRUE) %>%  ungroup() # 去掉基因的sum列重复行
#dist=distinct_rows %>% group_by(gene_name) %>% filter(n_distinct(sum) > 1) %>% ungroup() # 找出 sum 值不同的 gene_name
pan_core_info2= pan_core_info  %>% group_by(gene_name) %>% slice(1)
pan_core_info2$gene_name <- gsub("Msazm4", "Msa", pan_core_info2$gene_name)
combine_info=all_gene %>% left_join(pan_core_info2, by=c('gene'='gene_name')) 
##six group genes
combine_info$Type2=combine_info$Type

combine_info <- combine_info %>%
  mutate(Type2 = if_else(
    !is.na(sum) & Type == "four_copy" & sum == "Core", 
    "Fourfold_core", 
    Type2
  ))
combine_info <- combine_info %>%
  mutate(Type2 = if_else(
    !is.na(sum) & Type == "one_copy" & sum == "Dispensable", 
    "Onefold_essential", 
    Type2
  ))
group_info <- combine_info %>%
  group_by(Type2) %>%
  summarise(count = n(),
            proportion = count / nrow(combine_info)
  )
group_info

write.table(combine_info,"all_gene_pan_core_copy_number_gene_info_202473.txt",row.names = F,col.names = T,sep = "\t",quote = F)
combine_info=read.table("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info\\all_gene_pan_core_copy_number_gene_info_202473.txt",head = T,sep = "\t")

##new figure, Figure3B
####由于涉及到重复基因信息，因此将4,3,2份拷贝基因保留一份，然后计算不同基因的比例
unique_all_gene=all_gene
unique_all_gene$keep=NA
unique_all_gene$keep[unique_all_gene$gene %in% unique_gene_four$gene_name] <-"Yes"  ##unique four copy gene:22012
unique_all_gene$keep[unique_all_gene$gene %in% unique_gene_three$gene_name] <-"Yes" ##unique three copy gene:15264
unique_all_gene$keep[unique_all_gene$gene %in% unique_gene_two$gene_name] <-"Yes" ##unique two copy gene:23095
unique_all_gene$keep[unique_all_gene$gene %in% one_copy$gene] <-"Yes" ##unique one copy gene:49344
unique_all_gene=unique_all_gene[!is.na(unique_all_gene$keep),]
group_info <- unique_all_gene %>%
  group_by(Type) %>%
  summarise(count = n(),
            proportion = count / nrow(unique_all_gene)
  )
group_info
combine_info_fig=combine_info[combine_info$gene %in% unique_all_gene$gene,]


result <- combine_info_fig %>%
  group_by(Type, sum) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Type) %>%
  mutate(proportion = count / sum(count))
colnames(result)=c("Type","Group","count","proportion")
levels(result$Group) <- c(levels(result$Group), "Ungrouped")
result$Group[is.na(result$Group)] <- "Ungrouped"
library(stringr)
result <- result %>%
  mutate(Type = str_to_title(Type))
result$Type <- gsub("_", "-", result$Type)


result$Type <- factor(result$Type,levels = c("One-copy","Two-copy","Three-copy","Four-copy"),ordered = TRUE) 
result$Group <- factor(result$Group,levels = c("Ungrouped","Private","Dispensable","Softcore","Core"),ordered = TRUE)
p1=ggplot(data=result,aes(x = Type, y = count, group = Group, fill = Group)) +
  scale_fill_manual(breaks=c("Core","Softcore","Dispensable","Private","Ungrouped"),
                    values=c("#A4CB9E","#EDA1A4","#F9C08A","#9DB4CE","darkgray"))+  ##不同类型基因对应的颜色
  geom_bar(stat = "identity") +
  labs(x = "",y = "Number of genes") +
  theme(
    panel.background = element_blank(),
    axis.title=element_text(size=18,colour = 'black'), #坐标轴标题
    axis.text=element_text( size=18,colour = 'black'), #坐标轴标签
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size=0.5, colour = 'black'), #轴线
    plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.5, unit = "cm"))
ggsave(p1,filename = "bar_core_gene_copy_number_new_figure3b_unique_gene_new.pdf",width = 7,height = 6,dpi=300)
