###_________________________R CODE______________________________________________________
setwd("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info")
library(dplyr)
library(tidyr)
##get gene chr and pos information
pos_info=read.table("ZM4_V2.bed",header = F,sep = "\t")
colnames(pos_info)=c("Chr", "Start","end" , "gene","value","direction")

##1. extract gene paired info, example hap3 mapping to hap4
mydata=read.table("ZM4_V2_hap3_mapping_hap4.pep.blast",header = F,sep = "\t")
colnames(mydata)=c("qseqid", "sseqid",    "pident",    "length",    "mismatch",    "gapopen",    "qstart",    "qend",    "sstart",    "send",    "evalue",    "bitscore")
##filter evalue less than 1e-10
mydata=mydata[mydata$evalue<1e-10,]
#filter similarity >95%
mydata=mydata[mydata$pident>95,]
#filter coverage >95%
mydata["coverage"]=(1-mydata$mismatch/mydata$length)
mydata=mydata[mydata$coverage>0.95,]
mydata2 <- mydata %>% left_join(pos_info, by=c('qseqid'='gene'))  ##query gene ID
mydata3=mydata2 %>% left_join(pos_info, by=c('sseqid'='gene'))   ##match gene ID
mydata4 <- mydata3 %>% separate(Chr.x, into = c("Chr_quary", "hap_quary"), sep = "_")
mydata5 <- mydata4 %>% separate(Chr.y, into = c("Chr_match", "hap_match"), sep = "_")
mydata5 <- mydata5 %>% filter(Chr_quary == Chr_match)
mydata5=mydata5[order(mydata5$qseqid,mydata5$pident,decreasing = TRUE),]
result <- mydata5 %>%group_by(qseqid) %>% slice(1)
write.table(result,"ZM4_V2_hap3_mapping_hap4_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)

###2. extract four copied genes
##hap1比对到hap2-4之前已经比对过，没有新的基因对
##hap2比对到hap1,之前没有比对过，有新的基因对
##hap3比对到hap1-2,之前没有比对过，有新的基因对
##hap4之前没有比对过hap1-3,有新基因对
###合并第一次，第二次比对的结果
###combine hap1-4
hap12=read.table("ZM4_V2_hap1_mapping_hap2_gene.txt",header = T,sep = "\t")
hap13=read.table("ZM4_V2_hap1_mapping_hap3_gene.txt",header = T,sep = "\t")
hap14=read.table("ZM4_V2_hap1_mapping_hap4_gene.txt",header = T,sep = "\t")
hap123=hap12 %>% inner_join(hap13, by=c('qseqid'='qseqid'))
hap1234=hap123 %>% inner_join(hap14, by=c('qseqid'='qseqid'))
##remapped, hap4-1 four copy gene
hap41=read.table("ZM4_V2_remapped_hap4_mapping_hap1_gene.txt",header = T,sep = "\t")
hap42=read.table("ZM4_V2_remapped_hap4_mapping_hap2_gene.txt",header = T,sep = "\t")
hap43=read.table("ZM4_V2_remapped_hap4_mapping_hap3_gene.txt",header = T,sep = "\t")
hap412=hap41 %>% inner_join(hap42, by=c('qseqid'='qseqid'))
hap4123=hap412 %>% inner_join(hap43, by=c('qseqid'='qseqid'))
##combine
four_hap_gene=rbind(hap1234,hap4123)
write.table(four_hap_gene,"Four_haps_gene_paired_info.txt",row.names = F,col.names = T,sep = "\t",quote = F) ##export gene pair info

four_hap_gene$pair=1:nrow(four_hap_gene)
four_hap_gene_1=four_hap_gene[,c("qseqid","hap_quary.x","pair")]
four_hap_gene_2=four_hap_gene[,c("sseqid.x","hap_match.x","pair")]
four_hap_gene_3=four_hap_gene[,c("sseqid.y","hap_match.y","pair")]
four_hap_gene_4=four_hap_gene[,c("sseqid","hap_match","pair")]
##rename
colnames(four_hap_gene_1)=c("gene","hap","pair")
colnames(four_hap_gene_2)=c("gene","hap","pair")
colnames(four_hap_gene_3)=c("gene","hap","pair")
colnames(four_hap_gene_4)=c("gene","hap","pair")
combine_four_hap_gene1234=rbind(four_hap_gene_1,four_hap_gene_2,four_hap_gene_3,four_hap_gene_4)
##四拷贝成对信息
combine_four_hap_gene1234_sep=combine_four_hap_gene1234 %>% separate(gene, into = c("gene_name", "mRNA"), sep = ".t")
write.table(combine_four_hap_gene1234_sep,"combine_four_hap_gene1234_paired_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)

##remove duplicate info
unique_four_hap_gene1234 <- combine_four_hap_gene1234_sep %>%
  group_by(gene_name) %>%
  slice(which.min(pair)) %>%
  ungroup()  ##总79455个基因
##成对匹配
four_hap_paired_hap1=unique_four_hap_gene1234[unique_four_hap_gene1234$hap==1,] ##20845
four_hap_paired_hap2=unique_four_hap_gene1234[unique_four_hap_gene1234$hap==2,] ##19209
four_hap_paired_hap3=unique_four_hap_gene1234[unique_four_hap_gene1234$hap==3,] ##19271
four_hap_paired_hap4=unique_four_hap_gene1234[unique_four_hap_gene1234$hap==4,] ##20130
##唯一基因数量,22012
unique_gene_four <- unique_four_hap_gene1234 %>%
  group_by(pair) %>%
  slice(which.min(hap)) %>%
  ungroup()

##three copy, 124,134,123, 234
hap124=hap12 %>% inner_join(hap14, by=c('qseqid'='qseqid'))
hap134=hap13 %>% inner_join(hap14, by=c('qseqid'='qseqid'))
hap23=read.table("ZM4_V2_hap2_mapping_hap3_gene.txt",header = T,sep = "\t")
hap24=read.table("ZM4_V2_hap2_mapping_hap4_gene.txt",header = T,sep = "\t")
hap234=hap23 %>% inner_join(hap24, by=c('qseqid'='qseqid'))
##remapped
##three copy gene
hap423=hap42 %>% inner_join(hap43, by=c('qseqid'='qseqid'))
hap31=read.table("ZM4_V2_remapped_hap3_mapping_hap1_gene.txt",header = T,sep = "\t")
hap32=read.table("ZM4_V2_remapped_hap3_mapping_hap2_gene.txt",header = T,sep = "\t")
hap312=hap31 %>% inner_join(hap32, by=c('qseqid'='qseqid'))
three_hap=rbind(hap124,hap134,hap123,hap234,hap412,hap423,hap312) ##it may contains four copy genes, filter it in the following steps 
three_hap$pair=1:nrow(three_hap)
three_hap_1=three_hap[,c("qseqid","hap_quary.x","pair")]
three_hap_2=three_hap[,c("sseqid.x","hap_match.x","pair")]
three_hap_3=three_hap[,c("sseqid.y","hap_match.y","pair")]

##rename
colnames(three_hap_1)=c("gene","hap","pair")
colnames(three_hap_2)=c("gene","hap","pair")
colnames(three_hap_3)=c("gene","hap","pair")
combine_three_hap123=rbind(three_hap_1,three_hap_2,three_hap_3)
combine_three_hap123_gene_sep=combine_three_hap123 %>% separate(gene, into = c("gene_name", "mRNA"), sep = ".t")
combine_three_hap123_gene=combine_three_hap123_gene_sep[!combine_three_hap123_gene_sep$gene_name %in% combine_four_hap_gene1234_sep$gene_name,] #filter
three_hap_gene=three_hap[three_hap$pair %in% combine_three_hap123_gene$pair,]
write.table(three_hap_gene,"Three_haps_gene_paired_info.txt",row.names = F,col.names = T,sep = "\t",quote = F) ##export gene pair info

##三拷贝成对信息
##remove duplicate info
unique_three_hap_gene123 <- combine_three_hap123_gene %>%
  group_by(gene_name) %>%
  slice(which.min(pair)) %>%
  ungroup()  ##总36330个基因
##成对匹配
three_hap_paired_hap1=unique_three_hap_gene123[unique_three_hap_gene123$hap==1,] ##9117
three_hap_paired_hap2=unique_three_hap_gene123[unique_three_hap_gene123$hap==2,] ##9604
three_hap_paired_hap3=unique_three_hap_gene123[unique_three_hap_gene123$hap==3,] ##8790
three_hap_paired_hap4=unique_three_hap_gene123[unique_three_hap_gene123$hap==4,] ##8819
##唯一基因数量,15264
unique_gene_three <- unique_three_hap_gene123 %>%
  group_by(pair) %>%
  slice(which.min(hap)) %>%
  ungroup()

##two copy gene
hap34=read.table("ZM4_V2_hap3_mapping_hap4_gene.txt",header = T,sep = "\t")
hap21=read.table("ZM4_V2_remapped_hap2_mapping_hap1_gene.txt",header = T,sep = "\t")
two_hap=rbind(hap12,hap13,hap14,hap23,hap24,hap34,hap41,hap42,hap43,hap31,hap32,hap21)
two_hap$pair=1:nrow(two_hap)
two_hap_1=two_hap[,c("qseqid","hap_quary","pair")]
two_hap_2=two_hap[,c("sseqid","hap_match","pair")]


##rename
colnames(two_hap_1)=c("gene","hap","pair")
colnames(two_hap_2)=c("gene","hap","pair")
combine_two_hap12=rbind(two_hap_1,two_hap_2)
combine_two_hap12_gene_sep=combine_two_hap12 %>% separate(gene, into = c("gene_name", "mRNA"), sep = ".t")

combine_two_hap12_gene=combine_two_hap12_gene_sep[!combine_two_hap12_gene_sep$gene_name %in% combine_four_hap_gene1234_sep$gene_name,]
combine_two_hap12_gene=combine_two_hap12_gene[!combine_two_hap12_gene$gene %in% combine_three_hap123_gene$gene_name,]

two_hap_gene=two_hap[two_hap$pair %in% combine_two_hap12_gene$pair,]
write.table(two_hap_gene,"Two_haps_gene_paired_info.txt",row.names = F,col.names = T,sep = "\t",quote = F) ##export gene pair info

##二拷贝成对信息
##remove duplicate info
unique_two_hap_gene12 <- combine_two_hap12_gene %>%
  group_by(gene_name) %>%
  slice(which.min(pair)) %>%
  ungroup()  ##总37344个基因
##成对匹配
two_hap_paired_hap1=unique_two_hap_gene12[unique_two_hap_gene12$hap==1,] ##8624
two_hap_paired_hap2=unique_two_hap_gene12[unique_two_hap_gene12$hap==2,] ##8828
two_hap_paired_hap3=unique_two_hap_gene12[unique_two_hap_gene12$hap==3,] ##10580
two_hap_paired_hap4=unique_two_hap_gene12[unique_two_hap_gene12$hap==4,] ##9312
##二拷贝唯一基因数量,23095
unique_gene_two <- unique_two_hap_gene12 %>%
  group_by(pair) %>%
  slice(which.min(hap)) %>%
  ungroup()
##one copy gene check
all_gene=read.table("E:\\Genome\\Comparative_genomics\\RNA_seq_results\\top10_percentage_gene_info\\ZM4_V2.bed",header = T,sep = "\t")

##一拷贝基因
one_copy=all_gene[!all_gene$gene %in% combine_four_hap_gene1234_sep$gene_name, ] ##所有四拷贝总量：79455,39.2%
one_copy=one_copy[!one_copy$gene %in% combine_three_hap123_gene$gene_name, ] ##所有三拷贝总量：36330, 17.9%
one_copy=one_copy[!one_copy$gene %in% combine_two_hap12_gene$gene_name, ]  ##所有二拷贝基因总量：37344, 18.4%
##49344个一拷贝基因,24.4%
write.table(one_copy,"One_hap_gene_paired_info.txt",row.names = F,col.names = T,sep = "\t",quote = F) ##export gene pair info

##combine gene information
all_gene$Type=NA
all_gene$Type[all_gene$gene %in% combine_four_hap_gene1234_sep$gene_name] <-"four_copy"
all_gene$Type[all_gene$gene %in% combine_three_hap123_gene$gene_name] <-"three_copy"
all_gene$Type[all_gene$gene %in% combine_two_hap12_gene$gene_name] <-"two_copy"
all_gene$Type[all_gene$gene %in% one_copy$gene] <-"one_copy"
##计算比例
result <- all_gene %>%
  group_by(Type) %>%  # 按照 Type 列进行分组
  summarise(
    count = n(),  # 计算每一类的数量
    proportion = n() / nrow(all_gene)  # 计算每一类的比例
  )
result
write.table(all_gene,"ZM4_V2_gene_copy_info_202473_final_use.txt",row.names = F,col.names = T,sep = "\t",quote = F)
