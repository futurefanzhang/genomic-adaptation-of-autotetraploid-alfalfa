##1, blast each haplotype to other haplotypes
makeblastdb -dbtype prot -in  ZM4_V2_hap2_use.pep.fa -out ZM4_V2_hap2_database
blastp -query ZM4_V2_hap1_use.pep.fa -db ZM4_V2_hap2_database -out ZM4_V2_hap1_mapping_hap2.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 12
blastp -query ZM4_V2_hap1_use.pep.fa -db ZM4_V2_hap3_database -out ZM4_V2_hap1_mapping_hap3.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 12
blastp -query ZM4_V2_hap1_use.pep.fa -db ZM4_V2_hap4_database -out ZM4_V2_hap1_mapping_hap4.pep.blast  -outfmt 6  -num_alignments 5  -num_threads 12
……
##2, extract gene info(chr,pos)
python -m jcvi.formats.gff bed --type=mRNA --key=ID ZM4_V2.gff3 > ZM4_V2.bed ##MCScanX软件，提取蛋白对应的坐标
##3,using R do filter analysis
##gene_copy_number.R
