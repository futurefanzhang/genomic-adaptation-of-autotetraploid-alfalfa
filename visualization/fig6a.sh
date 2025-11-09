python -m jcvi.compara.synteny mcscan ZM4_V2_hap1.bed ZM4_V2_hap1.ZM4_V2_hap2.lifted.anchors --iter=1 -o ZM4_V2_hap1.ZM4_V2_hap2.i1.blocks ##更加细节的基因层面共线性
python -m jcvi.formats.base join ZM4_V2_hap1.Mara.i1.blocks ZM4_V2_hap1.ZM4_V2_hap2.i1.blocks ZM4_V2_hap1.Mca_long.i1.blocks --noheader | cut -f1,2,4,6 > seven_species1.blocks ##把不同的blocks合并起来，这块一次只能合并四个基因组信息，分两次进行
python -m jcvi.formats.base join  ZM4_V2_hap1.ZM4_V2_hap3.i1.blocks ZM4_V2_hap1.ZM4_V2_hap4.i1.blocks ZM4_V2_hap1.Mlup.i1.blocks  --noheader | cut -f1,2,4,6 > seven_species2.blocks
python -m jcvi.formats.base join ZM4_V2_hap1.Mt_5.i1.blocks ZM4_V2_hap1.Mpoly.i1.blocks ZM4_V2_hap1.Mru_zhiwusuo.i1.blocks --noheader | cut -f1,2,4,6 > seven_species3.blocks
##filter genes within up and down of Msa107765.t01
awk -v pattern="Msa107765.t01" -v n=30 '{lines[NR] = $0}
$0 ~ pattern {
    start = NR - n
    if (start < 1) start = 1
    end = NR + n
    for (i = start; i <= NR; i++) {
        if (i in lines) print lines[i]
    }
    for (i = 1; i <= n; i++) {
        if ((getline line) > 0) {
            print line
        } else {
            break
        }
    }
    exit
}' seven_species3.blocks > use_seven_species3.blocks

##block file
Msa107761.t01   Mara39886.t01   Msa113773.t01   Mca0397900      Msa119643.t01   Msa125464.t01   Mlup43583.t01   MtrunA17Chr5g0400411    mRNA.Mpo5G4710  evm.model.fragScaff_scaffold_66_pilon.18
Msa107762.t01   Mara39886.t01   Msa113774.t01   Mca0223640      Msa119644.t01   Msa125464.t01   Mlup43583.t01   MtrunA17Chr5g0400411    mRNA.Mpo5G4700  evm.model.fragScaff_scaffold_66_pilon.18
Msa107763.t01   .       .       .       .       .       .       .       .       .
Msa107763.t02   .       .       .       .       .       .       .       .       .
Msa107763.t04   .       .       .       .       .       .       .       .       .
Msa107764.t01   Mara39884.t01   Msa113775.t01   Mca0415500      Msa119645.t01   Msa125464.t01   Mlup43583.t01   MtrunA17Chr5g0400421    mRNA.Mpo5G4710  evm.model.fragScaff_scaffold_66_pilon.18
magenta*Msa107765.t01   Mara39883.t01   Msa113776.t01   Mca0415510      Msa119646.t01   Msa125466.t01   Mlup43586.t01   MtrunA17Chr5g0400431    mRNA.Mpo5G4720  evm.model.fragScaff_scaffold_66_pilon.17
Msa107766.t01   .       Msa113777.t01   .       .       .       .       .       .       .
Msa107767.t01   Mara39882.t01   Msa113778.t01   Mca0415520      Msa119647.t01   Msa125468.t01   Mlup43587.t01   MtrunA17Chr5g0400451    mRNA.Mpo5G4730  evm.model.fragScaff_scaffold_66_pilon.16
Msa107768.t01   Mara39881.t01   Msa113779.t01   Mca0415530      Msa119648.t01   Msa125469.t01   Mlup43588.t01   MtrunA17Chr5g0400461    mRNA.Mpo5G4740  evm.model.fragScaff_scaffold_66_pilon.15

##layout file
# x,   y, rotation,     ha,     va, color, ratio,            label
0.5, 0.65,        0, left,    center,      darkorange,     1, ZM4_V2_hap1 Chr5
0.2, 0.8,        30, center, top,      black,    .5, Mara Chr5
0.5, 0.55,        0, left,    center,      darkorange,     1, ZM4_V2_hap2 Chr5
0.5, 0.8,        0, center, top,      black,    .5, Mca_long Chr5
0.5, 0.45,        0, left,    center,      darkorange,     1, ZM4_V2_hap3 Chr5
0.5, 0.35,        0, left,    center,     darkorange,     1, ZM4_V2_hap4 Chr5
0.2, 0.2,        -30, center, bottom,      black,    .5, Mlup Chr5
0.5, 0.2,        0, center, bottom,      black,    .5, Mt_5 Chr5
0.8, 0.2,        30, center, bottom,      black,    .5, Mpoly Chr5
0.8, 0.8,        -30, center, top,      black,    .5, Mru_zhiwusuo Chr5
# edges
e, 0, 1
e, 0, 2
e, 0, 3
e, 0, 9
e, 0, 9
e, 2, 4
e, 4, 5
e, 5, 6
e, 5, 7
e, 5, 8

##Run the following command in the terminal:
python -m jcvi.formats.bed merge ZM4_V2_hap1.bed Mara.bed ZM4_V2_hap2.bed Mca_long.bed ZM4_V2_hap3.bed ZM4_V2_hap4.bed Mlup.bed Mt_5.bed Mpoly.bed Mru_zhiwusuo.bed -o use_seven_species_all.bed
python -m jcvi.graphics.synteny use_seven_species_all.blocks use_seven_species_all.bed use_seven_species_all.layout
