# genomic-adaptation-of-autotetraploid-alfalfa

This is the script used for genomic analysis of climate adaptation of alfalfa

### step1, genome assembly of alfalfa, including tetraploid variants calling

you can see the detail information from here: 
```
snp_and_indel_calling.sh
freebayes_variants_calling.sh

```
### step2, Identification of Gene Copy Number Variations

This part includes all related content of DNA and RNA-seq related analysis,
```
phylogenetic tree: phylogenetic_tree.sh, phylogenetic_tree_plot.R, vcf2other.py
PCA and IBD: PCA_and_IBD.sh, PCA.R, IBD.R, lecture06_07_add_id.pl
