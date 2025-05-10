#extract information in orthofinder
cut -f1,13 /data/home/zhangfan/compare_genomoe_medicago/genespace_new/orthofinder/Results_May09/Orthogroups/Orthogroups.tsv > ZM4_V2_gene_family_detail.txt ##提取目标基因组的信息
python change_gene_family.py ZM4_V2_gene_family_detail.txt ZM4_V2_gene_family_detail_change_format.txt ##将每行多基因拆分为每行一个基因的格式
##follow steps in R
setwd("E:\\Genome\\Comparative_genomics\\core_dispensable_genes_v2/")
##gene pair info
mydata=read.table("pan_gene_info.txt",row.names = 1,head = T,sep = "\t")
gene_detail=read.table("ZM4_V2_gene_family_detail_change_format.txt",head = T,sep = "\t")
mydata$Orthogroup=rownames(mydata)
mydata_use=subset(mydata, select = c(Orthogroup,sum))  ##select: Orthogroup,sum
library(dplyr)
all_gene_info <- gene_detail %>% left_join(mydata_use, by=c('Orthogroup'='Orthogroup'))
write.table(all_gene_info,"all_gene_orthogroup_pan_gene_info.txt",row.names = F,col.names = T,sep = "\t",quote = F)

counts <- all_gene_info %>%
  group_by(sum) %>%          # 按 gene 列分组
  summarize(count = n()) %>% # 计算每组的计数
  mutate(proportion = count / sum(count)) # 计算比例
counts 
##统计每一类基因的占比，我的结果核心基因占比35.8%， 私有基因10.8%
