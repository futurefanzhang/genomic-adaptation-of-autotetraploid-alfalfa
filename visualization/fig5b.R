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
  group_by(Type2, region) %>%
  summarize(
    mean_value = mean(mean_value, na.rm = TRUE))%>%
  mutate(
    proportion = mean_value / sum(mean_value)
  ) ##generate the mean deleterious SNP number of each chr
result1 
result <- result1 %>%
  group_by(region) %>%
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
ggsave(p1,filename = "Number of deleterious variance_core_gene_region2.pdf",width = 12,height = 6,dpi=300)
