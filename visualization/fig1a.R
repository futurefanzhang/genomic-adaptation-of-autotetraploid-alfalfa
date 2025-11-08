prepare data
conda activate Genome ##需要安装R,orthonfinder,jcvi
python -m jcvi.formats.gff bed --type=mRNA --key=ID A_M411.gff3.gz -o A_M411.bed
##analysis in R
setwd("/data/home/zhangfan/compare_genomoe_medicago/genespace1") ##修改分析文件的路径
library(GENESPACE)
gpar <- init_genespace(
  wd = "/data/home/zhangfan/compare_genomoe_medicago/genespace1", ##文件存在的路径
  path2mcscanx = "/data/home/zhangfan/genome/co_linear/MCScanX-master/")##检查数据格式，并创建分析对应的文件夹
gpar <- run_genespace(gsParam = gpar) ##前面运行成功的结果会跳过
##modify figure
load('/data/home/zhangfan/compare_genomoe_medicago/genespace1/results/gsParams.rda',verbose = TRUE) ##载入之前分析的结果
ggthemes <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))  ##修改背景的黑色为白色
ripd <- plot_riparian(
  gsParam,
  addThemes = ggthemes,chrBorderCol="black", ##载入白色背景，修改染色体边框为黑色
  refGenome = "ZM4_V2", ##设置参考基因组
genomeIDs =c("Mru_landa","Mru_zhiwusuo","Mlup","Mara","Mpoly","Mt_HM078","Mt_5","Mt_R108","Mca.landa","Mca_long","xinjiangdaye","zhongmu1","ZM4_V2"), ##修改不同基因组出现的顺序，从左到右代表图形上面从下到上的排列方式
pdfFile = "replot.pdf" ##导出的图形文件
)
