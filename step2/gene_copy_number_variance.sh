##1, blast each haplotype to other haplotypes
makeblastdb -dbtype prot -in  ZM4_V2_hap2_use.pep.fa -out ZM4_V2_hap2_database
blastp -query ZM4_V2_hap1_use.pep.fa -db ZM4_V2_hap2_database -out ZM4_V2_hap1_mapping_hap2.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 12
blastp -query ZM4_V2_hap1_use.pep.fa -db ZM4_V2_hap3_database -out ZM4_V2_hap1_mapping_hap3.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 12
blastp -query ZM4_V2_hap1_use.pep.fa -db ZM4_V2_hap4_database -out ZM4_V2_hap1_mapping_hap4.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 12
……
##2, extract gene info(chr,pos)
python -m jcvi.formats.gff bed --type=mRNA --key=ID ZM4_V2.gff3 > ZM4_V2.bed ##MCScanX软件，提取蛋白对应的坐标
##3,using R do filter
_________________________R CODE______________________________________________________
setwd("E:\\Genome\\Comparative_genomics\\ZM4_1-4copy_gene_info")
library(dplyr)
library(tidyr)
##get gene chr and pos information
pos_info=read.table("ZM4_V2.bed",header = F,sep = "\t")
colnames(pos_info)=c("Chr", "Start","end" , "gene","value","direction")

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

###combine hap1-4
hap12=read.table("ZM4_V2_hap1_mapping_hap2_gene.txt",header = T,sep = "\t")
hap13=read.table("ZM4_V2_hap1_mapping_hap3_gene.txt",header = T,sep = "\t")
hap14=read.table("ZM4_V2_hap1_mapping_hap4_gene.txt",header = T,sep = "\t")
hap123=hap12 %>% inner_join(hap13, by=c('qseqid'='qseqid'))
hap1234=hap123 %>% inner_join(hap14, by=c('qseqid'='qseqid'))
write.table(hap1234,"ZM4_V2_hap1-4mapping_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)

##three copy, 124,134,123, 234
hap124=hap12 %>% inner_join(hap14, by=c('qseqid'='qseqid'))
hap134=hap13 %>% inner_join(hap14, by=c('qseqid'='qseqid'))
hap23=read.table("ZM4_V2_hap2_mapping_hap3_gene.txt",header = T,sep = "\t")
hap24=read.table("ZM4_V2_hap2_mapping_hap4_gene.txt",header = T,sep = "\t")
hap234=hap23 %>% inner_join(hap24, by=c('qseqid'='qseqid'))
three_hap=rbind(hap123,hap124,hap134,hap234)
gene1=as.data.frame(unique(three_hap$qseqid))
colnames(gene1)="three_gene"
gene2=as.data.frame(unique(three_hap$sseqid.x))
colnames(gene2)="three_gene"
gene3=as.data.frame(unique(three_hap$sseqid.y))
colnames(gene3)="three_gene"
three_gene=rbind(gene1,gene2,gene3)
four_gene=read.table("four_copyhap1-4gene.txt",header = T,sep = "\t") ##根据hap1234的结果，处理成一列，84522, 总：247127,比例：34.20%
three_gene <- as.data.frame(three_gene[!three_gene$three_gene %in% four_gene$four_gene,])
colnames(three_gene)="three_gene" ##44419， 总：247127,比例：17.97%
write.table(three_gene,"ZM4_V2_three_gene_mapping_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)

##two copy gene
hap34=read.table("ZM4_V2_hap3_mapping_hap4_gene.txt",header = T,sep = "\t")
two_hap=rbind(hap12,hap13,hap14,hap23,hap24,hap34)
gene1=as.data.frame(unique(two_hap$qseqid))
colnames(gene1)="two_gene"
gene2=as.data.frame(unique(two_hap$sseqid))
colnames(gene2)="two_gene"
two_gene=rbind(gene1,gene2)
two_gene <- as.data.frame(two_gene[!two_gene$two_gene %in% four_gene$four_gene,])
two_gene <- as.data.frame(two_gene[!two_gene$two_gene %in% three_gene$three_gene,])
colnames(two_gene)="two_gene" ##45253， 总：247127,比例：18.31%
write.table(two_gene,"ZM4_V2_two_gene_mapping_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)
_________________________R CODE______________________________________________________
##step4, mapping rest genes
cat four_copyhap1-4gene.txt ZM4_V2_three_gene_mapping_gene.txt ZM4_V2_four_gene_mapping_gene.txt > ZM4_V2_234_copy_gene.txt
seqkit grep -v -f ZM4_V2_234_copy_gene.txt ZM4_V2.pep.fa > ZM4_V2_one_copy.pep.fa
blastp -query ZM4_V2_one_copy.pep.fa -db ZM4_V2_four_copy_database -out ZM4_V2_one_copy_mapping_four_hap.pep.blast  -outfmt 6  -num_alignments 10  -num_threads 30
##step5, double check copy number, Rcode
_________________________R CODE______________________________________________________
##one copy gene check
hap_one=read.table("ZM4_V2_one_copy_mapping_four_hap.pep.blast",header = T,sep = "\t")
colnames(hap_one)=c("qseqid", "sseqid",    "pident",    "length",    "mismatch",    "gapopen",    "qstart",    "qend",    "sstart",    "send",    "evalue",    "bitscore")
##filter evalue less than 1e-10
hap_one=hap_one[hap_one$evalue<1e-10,]
#filter similarity >90%
hap_one=hap_one[hap_one$pident>95,]
#filter coverage >90%
hap_one["coverage"]=(1-hap_one$mismatch/hap_one$length)
hap_one=hap_one[hap_one$coverage>0.95,]
mydata2 <- hap_one %>% left_join(pos_info, by=c('qseqid'='gene'))  ##query gene ID
mydata3=mydata2 %>% left_join(pos_info, by=c('sseqid'='gene'))   ##match gene ID
mydata4 <- mydata3 %>% separate(Chr.x, into = c("Chr_quary", "hap_quary"), sep = "_")
mydata5 <- mydata4 %>% separate(Chr.y, into = c("Chr_match", "hap_match"), sep = "_")
mydata5 <- mydata5 %>% filter(Chr_quary == Chr_match)
mydata5=mydata5[order(mydata5$qseqid,mydata5$pident,decreasing = TRUE),]
mydata5 <- mydata5 %>% filter(hap_quary != hap_match)
result <- mydata5 %>%group_by(qseqid) %>% slice(1)
write.table(result,"ZM4_V2_one_copy_mapping_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)
_________________________R CODE______________________________________________________
##step6, non mapped genes were consider as single copy gene
cat ZM4_V2_234_copy_gene.txt  ZM4_V2_one_copy_remapped_gene.txt |sed -e "s/\r//g" > ZM4_V2_more_than_onecpopy_gene.txt #检查是否有特殊字符^M，并去掉，注意：这部分可能有重复基因
grep ">" ZM4_V2.pep.fa | sed "s/>//g" | sed -e "s/\r//g" > ZM4_V2_all_gene.txt
awk 'FNR==NR{a[$0];next} !($0 in a)' ZM4_V2_more_than_onecpopy_gene.txt ZM4_V2_all_gene.txt > ZM4_V2_single_copy_gene.txt ##保留的单拷贝基因，比例为25.69%
##step7,remapp new mapped genes
seqkit grep -f remapped_hap1_gene.txt ZM4_V2.pep.fa >  ZM4_V2_remapped_hap1_gene.pep.fa
seqkit grep -f remapped_hap2_gene.txt ZM4_V2.pep.fa >  ZM4_V2_remapped_hap2_gene.pep.fa
seqkit grep -f remapped_hap3_gene.txt ZM4_V2.pep.fa >  ZM4_V2_remapped_hap3_gene.pep.fa
seqkit grep -f remapped_hap4_gene.txt ZM4_V2.pep.fa >  ZM4_V2_remapped_hap4_gene.pep.fa
blastp -query ZM4_V2_remapped_hap3_gene.pep.fa -db ZM4_V2_hap2_database -out ZM4_V2_remapped_hap3_mapping_hap2.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 10
blastp -query ZM4_V2_remapped_hap3_gene.pep.fa -db ZM4_V2_hap4_database -out ZM4_V2_remapped_hap3_mapping_hap4.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 10
……
##step8, refilter in R, same as step3
