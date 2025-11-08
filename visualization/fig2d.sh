#prepare density file
perl /data/home/zhangfan/compare_genomoe_medicago/PGGB/pggb/extract_SV_from_vcf.pl 0.vcf chr1_sv.vcf #ref或者alt超过50bp就认为是SV
perl /data/home/zhangfan/compare_genomoe_medicago/PGGB/pggb/calculate_SV_density_from_vcf.pl chr1_sv.vcf chr1_sv_density_1Mb.txt 1000000 ##输出文件为chr1_sv_density_1Mb.txt，区间大小设置为1Mb

##figures, using conda software circosplot
conda create -n circosplot circos
micromamba activate circosplot
cd /data/home/zhangfan/compare_genomoe_medicago/circosplot/circos_fig
把上面的四个统计结果放进一个data文件夹里面（如下图），再准备一个基因组信息文件karyotype.txt
circos -conf main_use.conf ##运行程序只有这一步，但是需要把下面配置文件都准备好
修改main_use.conf文件信息
karyotype=data/karyotype.txt    #基因组信息文件（如下图），文件放在data文件夹下
#基因组信息文件格式：
#chr - ID LABEL START END COLOR

chromosomes_units=1000000          #设置单位，1M，后面1u代表1M
chromosomes_display_default=yes    #展示所有染色体，如果为no，则需要指定chromosome参数

<<include conf/1.demoideogram.conf>>      #导入染色体配置参数，在conf文件夹下
<<include conf/2.demoticks.conf>>         #导入刻度配置参数，在conf文件夹下

<plots>                           #scatter line hist都属于plots, 放到一个<plots>中
<<include conf/heatmap.conf>>      #导入密度图信息，这个里面包括SNP,SV等四类变异信息，在conf文件夹下设置

</plots>

<image>
angle_offset*=-82 ##圈图开口大小
<<include  etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>

##配置文件1.demoideogram.conf如下，主要修改染色体条数就好
<ideogram>show=yes    #是否显示<spacing>default=0.005r   #设置圈图中染色体之间的空隙大小，以下设置为每个空隙大小为周长的 0.5%<pairwise chr1 chr8>spacing=5r</pairwise></spacing>radius=0.9r     #设定半径，染色体条在离圆心的 90% 处thickness=40p   #染色体条的宽度，可以使用 r（比例关系） 或 p（像素）作为单位fill=yes        #是否填充颜色。填充的颜色取决于 karyotype 指定的文件的最后一列#设定轮廓的颜色及其厚度。如果没有设置参数或设定其厚度为0，则没有轮廓stroke_color     = dgreystroke_thickness = 2p#设定是否显示 label，对应karyotype文件第4列。如果其值为 yes，则必须要有label_radius 参数来设定位置，否则会报错show_label=yes     label_radius=1.07r #label的半径label_font=default #label的字体label_size=50      #label的字体大小label_parallel=yes #label的字体显示方向，yes易于浏览show_bands=yes    #展示染色体组型fill_bands=yes    #染色体组型填充</ideogram>
配置文件2.demoticks.conf如下，基本无需修改

show_ticks=yes         #是否显示刻度线
show_tick_labels=yes   #是否显示刻度线label
<ticks>
radius=1.0r   #刻度线的半径
color=black    #刻度线的颜色
thickness=2p   #刻度线的宽度
multiplier=1e-6 #对显示的labels进行处理，对坐标乘以1e-6再显示，表示1M
format=%d       #整数格式,f for float
<tick>
spacing=1u     # 小刻度，每1个单位标注
size=10p        #刻度长度为5p
</tick>
<tick>
spacing=5u     #大刻度，每5个单位标注
size=20p       #刻度长度为10p
show_label=yes #显示刻度label
label_size=30p #label的大小
label_offset=10p #显示label的位置，一般与刻度长度一样就可以
format=%d      #显示格式为整数
</tick>
</ticks>

配置文件heatmap.conf设置

<plot>
type    = heatmap
# default file for all tracks
file             = data/use_snp_information.txt.binCount ##最外圈为SNP
r1=0.99r     #绘图外径
r0=0.90r     #绘图内径
# a 9 color diverging spectral palette specified using a color list name
color  = blues-7-seq-1,blues-7-seq-2,blues-7-seq-3,blues-7-seq-4,blues-7-seq-5  ##设置颜色
scale_log_base = 1.5
stroke_thickness = 0.001p
stroke_color     = white
min              = 0
max              = 1600 ##设置bincount的最大值，这个需要根据文件里面的信息确定一下
</plot>


<plot>
type    = heatmap
# default file for all tracks
file             = data/use_indel_information.txt.binCount ##第二圈为indel
r1=0.89r     #绘图外径
r0=0.80r     #绘图内径
# a 9 color diverging spectral palette specified using a color list name
color  = vlyellow,lyellow,greens-7-seq-4,greens-7-seq-5,greens-7-seq-6 ##设置颜色信息
scale_log_base = 1
stroke_thickness = 0.001p
stroke_color     = white
min              = 0
max              = 360  ##设置bincount的最大值，这个需要根据文件里面的信息确定一下
</plot>


<plot>
type    = heatmap
# default file for all tracks
file= data/use_sv_information.txt.binCount  ##第三圈为SV信息
r1=0.79r     #绘图外径
r0=0.70r     #绘图内径
# a 9 color diverging spectral palette specified using a color list name
color  = vvlyellow,vlyellow,gnbu-7-seq-3,gnbu-7-seq-4,gnbu-7-seq-5
scale_log_base = 1
stroke_thickness = 0.001p
stroke_color     = white
min              = 0
max              = 35 ##设置bincount的最大值，这个需要根据文件里面的信息确定一下
</plot>


<plot>
type    = heatmap
# default file for all tracks
file= data/SNP_freebayes_info.txt.binCount  ##第四圈为freebayes对应的SNP和indel变异信息
r1=0.69r     #绘图外径
r0=0.60r     #绘图内径
# a 9 color diverging spectral palette specified using a color list name
color  = vvlblue,vlblue,pubu-9-seq-5,pubu-9-seq-6,pubu-9-seq-7
scale_log_base = 1
stroke_thickness = 0
stroke_color     = white
min              = 0
max              = 3600 ##设置bincount的最大值，这个需要根据文件里面的信息确定一下
</plot>
