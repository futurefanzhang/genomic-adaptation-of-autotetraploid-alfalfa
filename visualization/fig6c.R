#https://dbsloan.github.io/TS2019/exercises/r_figure_drawing.html
#colinear analysis
setwd("E:\\Genome\\co_linear\\region_gene_structure")
genes=read.delim("genes.txt")

species1_min = min(genes[which(genes$Species==1), "Start"])
species1_max = max(genes[which(genes$Species==1), "End"])
species1_length = species1_max - species1_min + 1

species2_min = min(genes[which(genes$Species==2), "Start"])
species2_max = max(genes[which(genes$Species==2), "End"])
species2_length = species2_max - species2_min + 1

species3_min = min(genes[which(genes$Species==3), "Start"])
species3_max = max(genes[which(genes$Species==3), "End"])
species3_length = species3_max - species3_min + 1

species4_min = min(genes[which(genes$Species==4), "Start"])
species4_max = max(genes[which(genes$Species==4), "End"])
species4_length = species4_max - species4_min + 1

largestLength = max (species1_length, species2_length,species3_length,species4_length)

boxHeight = 0.04
height4 = 0.2
height3 = 0.4
height2 = 0.6
height1 = 0.8
pdf("regional_gene_structure_bin_30bp.pdf", width = 8, height = 4) 
plot(1, type = "n", xlab = "", ylab = "",  xaxt = "n",yaxt = "n", bty = "n",xlim = c(-200, 2000),ylim = c(0, 1))
axis(1, at = seq(0, 2000, by = 500),col.axis = "black", lwd = 2)

segments(1, height1, species1_length, height1, lwd=3)
text(-150, height1, labels = "Hap1",adj=0.5)
segments(1, height2, species2_length, height2, lwd=3)
text(-150, height2, labels = "Hap2",adj=0.5)
segments(1, height3, species3_length, height3, lwd=3)
text(-150, height3, labels = "Hap3",adj=0.5)
segments(1, height4, species4_length, height4, lwd=3)
text(-150, height4, labels = "Hap4",adj=0.5)

connections=read.delim("connections.txt")
for (i in 1:dim(connections)[1]){
  left1 = (connections[i, "StartA"] )
  right1 = (connections[i, "EndA"] )
  left2 = (connections[i, "StartB"] )
  right2 = (connections[i,"EndB"] )
  left3 = (connections[i, "StartC"] )
  right3 = (connections[i,"EndC"] )
  left4 = (connections[i, "StartD"] )
  right4 = (connections[i,"EndD"] )
  polygon(c(left1, left2, right2, right1), c(height1, height2, height2, height1), col = "gray85", border = NA)
  polygon(c(left2, left3, right3, right2), c(height2, height3, height3, height2), col = "gray85", border = NA)
  polygon(c(left3, left4, right4, right3), c(height3, height4, height4, height3), col = "gray85", border = NA)
  
} 

##Gene structure
genes=genes[genes$type!="gene" & genes$type!="intron",]
for (i in 1:dim(genes)[1]){
  #Define the "left" and "right" x-coordinates.
  if (genes[i, "Species"] == 1){
    left = (genes[i,"Start"])
    right = (genes[i,"End"])
    height = height1
  }else if (genes[i, "Species"] == 2){
    left = (genes[i,"Start"])
    right = (genes[i,"End"])
    height = height2
  }else if (genes[i, "Species"] == 3){
    left = (genes[i,"Start"])
    right = (genes[i,"End"])
    height = height3
  }else if (genes[i, "Species"] == 4){
    left = (genes[i,"Start"])
    right = (genes[i,"End"])
    height = height4
  }
  
  #we have already defined the gene colors in our input data file.
  geneFillColor=as.character(genes[i,"Color"]);
    polygon (c(left,  left,right, right), c(height-boxHeight/2, height+boxHeight/2,height+boxHeight/2,height-boxHeight/2), col = geneFillColor, lwd = 0.5)
}

##line
linear=read.delim("lines_bin_30bp.txt")

for (i in 1:dim(linear)[1]){
  #Define the "left" and "right" x-coordinates.
  if (linear[i, "Species"] == 1){
    left = (linear[i,"Pos"])
    right = (linear[i,"Pos"])
    height = height1
  }else if (linear[i, "Species"] == 2){
    left = (linear[i,"Pos"])
    right = (linear[i,"Pos"])
    height = height2
  }else if (linear[i, "Species"] == 3){
    left = (linear[i,"Pos"])
    right = (linear[i,"Pos"])
    height = height3
  }else if (linear[i, "Species"] == 4){
    left = (linear[i,"Pos"])
    right = (linear[i,"Pos"])
    height = height4
  }
  
  #we have already defined the gene colors in our input data file.
  geneFillColor=as.character(linear[i,"Color"]);
  
  #Draw arrows, distinguishing between forward and reverse oriented genes
    polygon (c(left, left, right, right), c(height-boxHeight/2,  height+boxHeight/2, height+boxHeight/2,height-boxHeight/2), border = "red", lwd = 0.2)
}
dev.off() 

