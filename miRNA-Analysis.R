#Reading cancer data from The Cancer Genome Atlas (TCGA) to variable "c" based on the user input. Specifically, each number is associated with a cancer (1 = Head and Neck Cancer, 2 = Kidney Cancer, 3 = Liver Cancer, 4 = Lung Cancer, 5 = Thyroid Cancer, and 6 = Breast Cancer).
type <- scan (what=numeric())
if (type==1) {
  c = read.table("file:///C:/Users/Siya/Downloads/HNSC.txt", header = T, row.names = 1,
                 stringsAsFactors = F)
} else if (type==2) {
  c = read.table("file:///C:/Users/Siya/Downloads/KIRC.txt", header = T, row.names = 1,
                 stringsAsFactors = F)
} else if (type==3) {
  c = read.table("file:///C:/Users/Siya/Downloads/LIHC.txt", header = T, row.names = 1,
                 stringsAsFactors = F)
} else if (type==4) {
  c = read.table("file:///C:/Users/Siya/Downloads/LUAD.txt", header = T, row.names = 1,
                 stringsAsFactors = F)
} else if (type==5) {
  c = read.table("file:///C:/Users/Siya/Downloads/THCA.txt", header = T, row.names = 1,
                 stringsAsFactors = F)
} else if (type==6) {
  c = read.table("file:///C:/Users/Siya/Downloads/BRCA.txt", header = T, row.names = 1,
                 stringsAsFactors = F)
}

#"c" consists of the expression of miRNAs in both controlled and cancer patients. Makes "c" a matrix so that it can be manipulated.
c
c = as.matrix(c)
is.matrix(c)

#Sorts data from "c" into "all_tumors" and "all_controls" dependingly if it is a patient that has cancer or not.
all_tumors = grep("Tumor",colnames(c),ignore.case = T)
all_tumors
all_controls = grep("Control",colnames(c),ignore.case = T)
all_controls

#Doctors can add their dataset into variable "f"
f = read.table("file:///C:/siya.txt", header = T, row.names = 1, stringsAsFactors = F)

#Combines TCGA dataset and doctor's dataset into variable "d"
d = cbind(c,f)

#Makes "d" a matrix so that it can be manipulated.
d = as.matrix(d)
is.matrix(d)

#Sorts data from "d" into "all_tumors" and "all_controls" dependingly if it is a patient that has cancer or not.
all_tumors = grep("Tumor",colnames(d),ignore.case = T)
all_controls = grep("Control",colnames(d),ignore.case = T)

#Heatmap that compares expression of miRNAs in controlled vs cancer samples so that the doctor can see the most differentially expressed miRNAs in their dataset
library(pheatmap)
d = d[rowSums(d) > 0,]
col.features = data.frame(Type = factor(rep(c("Tumor","Control"), each=length(all_tumors))))
rownames(col.features) = colnames(d)
pheatmap(log2(d + 1),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("darkblue", "white", "red2")) (256),
         border_color = NA,
         show_rownames = F,
         show_colnames = F,
         xlab="Sample",
         main="Heatmap",
         gaps_col=length(all_tumors),
         annotation_col=col.features,
         annotation_colors=list(Type = c(Tumor="red",Control="green")))

#Calculates the mean, standard deviation, and p-values (between controlled and cancer samples) of each miRNA in "d"
means = vector()
sds = vector()
pvals = vector()
for(row in rownames(d)) {
  means = c(means, mean(d[row,]))
  sds = c(sds, sd(d[row,]))
  pvals = c(pvals, t.test((d[row,all_tumors]), (d[row,all_controls]))$p.value)
}
stats = cbind(means, sds, pvals)
rownames(stats) = rownames(d)
colnames(stats) = c("Mean", "Standard Deviation", "P-Value")

#miRNAs ordered based on p-values from least to greatest to find the miRNAs that are most differentially expressed (between controlled and cancer samples) in "d"
ordered.pvals = order(pvals)
ordered.pvals
sorted.pvals = pvals[ordered.pvals]
sorted.pvals

#Finds the 30 most differentially expressed miRNAs in "d" and saves them to variable d.top30
rows.top30 = rownames(d)[ordered.pvals[1:30]]
rows.top30
d.top30 = d[rows.top30,]
d.top30 = d.top30[rowSums(d.top30) > 0,]
d.top30

#Finds the most differentially expressed miRNAs (pval < 0.05) in "d" and saves them to variable d.sig
d.sig = d[which(pvals < 0.05),]
d.sig

#Scatterplot is constructed to numerically and visually see the correlation between two miRNAs
plot(d["hsa-let-7a-2-3p",], d["hsa-miR-429",],
     col=c(rep("red", length(all_tumors)), rep("blue", length(all_tumors))),
     main="Sample Expression for hsa-let-7a-2-3p and hsa-miR-429",
     xlab="hsa-let-7a-2-3p",
     ylab="hsa-miR-429")
abline(lm(d["hsa-miR-429",]~d["hsa-let-7a-2-3p",]))
cor(d["hsa-let-7a-2-3p",], d["hsa-miR-429",])

#Heatmap that compares the expression of the 30 most differentially expressed miRNAs in controlled vs cancer samples
col.features = data.frame(Type = factor(rep(c("Tumor","Control"), each=length(all_tumors))))
rownames(col.features) = colnames(d.top30)
pheatmap(log2(d.top30 + 1),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("darkblue", "white", "red2")) (256),
         border_color = NA,
         show_rownames = F,
         show_colnames = F,
         main="Heatmap of 30 significant MiRNAs",
         gaps_col=length(all_tumors),
         annotation_col=col.features,
         annotation_colors=list(Type = c(Tumor="red",Control="green")))

#correlation coefficient calculated between the top 30 significant miRNAs and other miRNAs
corAll = cor(d[rows.top30,],d[31:60,])
corAll

#correlation coefficient calculated between the tumors and controls for the top 30 significant miRNAs
corTC = cor(d[rows.top30,all_tumors],d[rows.top30,all_controls])
corTC

#correlation coefficient calculated between the tumors and all the samples for the top 30 significant miRNAs
corAT = cor(d[rows.top30,],d[rows.top30,all_tumors])
corAT

#correlation coefficient calculated between the controls and all the samples for the top 30 significant miRNAs
corAC = cor(d[rows.top30,],d[rows.top30,all_controls])
corAC

#separates rows in "d.top30" based on whether their calculated ratio of means (d.ratio) is greater than 1 or less than or equal to 1. The rows with ratios greater than 1 are stored in "d.up", and the rows with ratios less than or equal to 1 are stored in "d.down".
d.up = data.frame()
d.down = data.frame()
for(i in 1:dim(d.top30)[1]) {
  d.tum = mean(d.top30[i,all_tumors], na.rm=T)
  d.con = mean(d.top30[i,all_controls], na.rm=T)
  d.ratio = d.tum / d.con
  if(!is.nan(d.ratio)) {
    if(d.ratio > 1) {
      d.up = rbind(d.up, d.top30[i,])
      rownames(d.up)[dim(d.up)[1]] = rownames(d.top30)[i]
    } else {
      d.down = rbind(d.down, d.top30[i,])
      rownames(d.down)[dim(d.down)[1]] = rownames(d.top30)[i]
    }
  }
}
colnames(d.up) = colnames(d.top30)
colnames(d.down) = colnames(d.top30)
d.up
d.down

#Reads in and displays file that contains a variety of miRNAs and the gene it is associated with
target = read.table("file:///C:/Users/Shipra/Desktop/Rcamp/h_mirtb_strongevid_target.txt", stringsAsFactors = F)
target = as.matrix(target)
target

#Allows doctor to input a miRNA and receive information about the gene the miRNA is associated with
MiRNA <- scan (what=numeric())
target[MiRNA,]

#Shows which genes are upregulated in "d"
target.up = vector()
for(rna in rownames(d.up)) {
  for(row in which(target[,1] == rna, arr.ind=T)) {
    target.up = c(target.up, target[row,2])
  }
}
target.up

#Shows which genes are downregulated in "d"
target.down = vector()
for(rna in rownames(d.down)) {
  for(row in which(target[,1] == rna, arr.ind=T)) {
    target.down = c(target.down, target[row,2])
  }
}
target.down
is.matrix(target.down)
target.down = as.matrix(target.down)
dim(target.down)
