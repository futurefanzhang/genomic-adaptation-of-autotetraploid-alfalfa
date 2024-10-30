##hifiasm
hifiasm -o primary_alfalfa -t80 --ul $ont $hifi
##juicer
bwa index /data/M1.whb/1new/juicer/reference/primary_alfalfap_ctg.fa

python ./juicer/misc/generate_site_positions.py DpnII alfalfa /juicer/reference/primary_alfalfap_ctg.fa


sh juicer.sh -t 20 -s DpnII -g primary_alfalfap_ctg.fa \
-d /data/juicer \
-D /software/juicer \
-z /juicer/reference/primary_alfalfap_ctg.fa \
-p /juicer/DpnII.chrom.sizes \
-y /juicer/XXhap1hap2_DpnII.txt

##3d-dna
run-assembly-visualizer.sh -q 0 -r M1.assembly merged_nodups.txt
##juicebox
##manually modification
