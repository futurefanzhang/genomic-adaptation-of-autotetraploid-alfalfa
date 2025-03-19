##This script contained two parts of the assembly, ZM4_V1.5 and ZM4_V2.0
##Part 1: ZM4_V1.5

##hifiasm, generate contig
hifiasm -o primary_alfalfa -t80 --ul $ont $hifi
##juicer
bwa index /data/M1.whb/1new/juicer/reference/primary_alfalfap_ctg.fa
python ./juicer/misc/generate_site_positions.py DpnII alfalfa /juicer/reference/primary_alfalfap_ctg.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' alfalfa_DpnII.txt > alfalfa_DpnII.chrom.sizes
sh juicer.sh -t 20 -s DpnII -g primary_alfalfap_ctg.fa \
-d /data/juicer \
-D /software/juicer \
-z /juicer/reference/primary_alfalfap_ctg.fa \
-p /juicer/alfalfa_DpnII.chrom.sizes \
-y /juicer/alfalfa_DpnII.txt
##ragtag

##3d-dna
run-assembly-visualizer.sh -q 0 -r M1.assembly merged_nodups.txt
##juicebox
##manually modification
##re3d-dna
bash /public/home/ac685n40ot/tools/3d-dna/run-asm-pipeline-post-review.sh -q 0 \
 -r alfalfa.review.assembly \
./reference/primary_alfalfap_ctg.fa \
./aligned/merged_nodups.txt > 3d.log
##after this step, we generated ZM4_V1.5 genome, including the alfalfa haplotype resolved 32 chromosomes (rename32.fa), and merged haplotype 8 chromosome(single8.fa)

##Part 2: ZM4_V2.0
##after we generate ZM4_V1.5 genome, the genome size is 2.5Gb, it means that part of genome haven't been assembled. so we used ZM4_V1.0 (https://figshare.com/s/fb4ba8e0b871007a9e6c) and ZM4_V1.5 to improve our genome.


