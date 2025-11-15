# genomic-adaptation-of-autotetraploid-alfalfa

This is the script used for genomic analysis of climate adaptation of alfalfa

### step1, genome assembly of alfalfa
you can see the detail information from related code. The content was corresponding to file name.
```
genome_assembly.sh
genomic_annotation.sh
```
### step2, pangenome construction and core gene identification.

This part includes pangenome construction related analysis.
```
pangenome.sh
orthofinder.R
core_gene_in_ZM4V2.sh
```

### step3, allelic copy number variation and six gene groups identification. 

This part includes gene copy number variation and six groups identification 
```
blast_four_haps.sh
core_gene_copy_number_overlap_group.R
gene_copy_number.R
```

### step4, GERP score. 

This part includes gerp score calculation and deleterious variants numbers in genes. 
```
gerp_score.sh
fourfold_core_genes_dsnp.R
deleterious_variants_gene_region.R
```

### step5, RNA-seq related analysis. 

This part includes FPKM, differentially expressed genes, top expressed genes and climate adaption related genes analysis. 
```
RNA_seq_analysis_fpkm_and_deg_gene.R
FPKM_value_gene_group.R
DEG_gene_group_proportion.R
Top_percent_FPKM.R
all_expressed_gene_FPKM.R
climate_gene_enrichment.R
top10_percent_gene_exp_detail_info.R
```

### step6, visualization. 

This part includes the script used for generate main figures and supplementary figures. 
```
fig1a.R
fig1b.R
....
```
