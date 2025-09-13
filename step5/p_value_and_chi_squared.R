##1.four allelic copies vs core genes, suppl. fig3 (sheet:pan_gene_copy_number)
Tetra_copy_core=11735
all_Tetra_copy_gene=22012
core_all=29191
non_redundant_genes=109715 
data_vector <- c(Tetra_copy_core, all_Tetra_copy_gene-Tetra_copy_core, core_all, non_redundant_genes-core_all)
data_table <- matrix(data_vector, nrow = 2, byrow = T)
rownames(data_table) <- c("Group_four_copies", "Group_all")
colnames(data_table) <- c("Core", "Non-Core")
chisq.test(data_table)
#results:X-squared = 6104.2, df = 1, p-value < 2.2e-16

##2. one allelic copies vs dispensable genes, suppl. fig3 (sheet:pan_gene_copy_number)
One_copy_dispensable=20562
all_one_copy_gene=49344
dispensable_all=37242
non_redundant_genes=109715
data_vector <- c(One_copy_dispensable, all_one_copy_gene-One_copy_dispensable, dispensable_all, non_redundant_genes-dispensable_all)
data_table <- matrix(data_vector, nrow = 2, byrow = T)
rownames(data_table) <- c("Group_one_copy", "Group_all")
colnames(data_table) <- c("Dispensable", "Non-Dispensable")
chisq.test(data_table)
#results:X-squared = 877.96, df = 1, p-value < 2.2e-16

##3. adaptation-associated genes enrichment
Tetra_copy_core_adaptation_associated_genes=836
all_adaptation_associated_genes=2429
Tetra_copy_core_genes=43582
all_genes=202473
data_vector <- c(Tetra_copy_core_adaptation_associated_genes, all_adaptation_associated_genes-Tetra_copy_core_adaptation_associated_genes, Tetra_copy_core_genes, all_genes-Tetra_copy_core_genes)
data_table <- matrix(data_vector, nrow = 2, byrow = T)
rownames(data_table) <- c("Group_adaptation", "Group_all")
colnames(data_table) <- c("Tetra_copy_core_genes", "Non-Tetra_copy_core_genes")
chisq.test(data_table)
##results:X-squared = 234.22, df = 1, p-value < 2.2e-16

##4.DEGs enrichment
##compare a group of number with one number
mydata=read.table("supply_table7_all_studies_DEG_proportion_of_gene_group.txt",head = T,sep = "\t")
mydata=mydata[mydata$Type2.group.=="Tetra_copy_core",]
hist(mydata$proportion) ##two peak, use median is more reasonable.
median(mydata$proportion) ##compare median differ with 21.5%, wilcox.test is suitable
wilcox.test(mydata$proportion, mu = 0.215) #genomic baseline proportion：21.5%
##results:V = 505, p-value = 0.0001934
