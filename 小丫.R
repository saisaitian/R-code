suppressMessages(require(WGCNA))  
suppressMessages(require(caret))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()###加快运行速度
enableWGCNAThreads(nThreads = 7)##设置电脑可用核数
options(stringsAsFactors = F)
data= read.csv("G:/小丫画图/mmc2.csv",header = T)
data=data[1:21]

dim(data)
fix(data)

datExpr0 = as.data.frame(t(data))
names(datExpr0) = data$Annotation.Divergence
datExpr0<-datExpr0[-1,]
rownames(datExpr0) = names(data)[-1]
fix(datExpr0)


##筛选方差前25%的基因##
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
dim(expro.upper)
datExpr0=expro.upper
str(datExpr0)
datExpr1<-data.matrix(datExpr0)
fix(datExpr1)


##这一步是为了减少运算量，因为一个测序数据可能会有好几万个探针，
##而可能其中很多基因在各个样本中的表达情况并没有什么太大变化，为了减少运算量，
##这里我们筛选方差前25%的基因。

gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK

### 样本数聚类
sampleTree = hclust(dist(datExpr1), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


datExpr = as.data.frame( datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


### 4 网络构建

net = blockwiseModules(datExpr, power = 16,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
## 一般 我们需要改一下power参数即可，其他可以选择默认

table(net$colors)

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]


###########
color<-unique(moduleColors)
color[1]
for (i  in 1:15) {
  y=assign(paste(color[i],"expr",sep = "."),datExpr[moduleColors==color[i]])
  write.table(y,paste(color[i],"xls",sep = "."),sep = "\t")
}













## 5 模块保守性分析
#模块保守性分析是用来检验参考网络中存在的某基因模块是否也存在于测试网络中，Z大于10，代表保守性好，大于2小于10
#代表保守性差，小于2代表无保守型，即参考网络中的模块在测试网络中并不存在，该模块是其特有的



inTraining <- createDataPartition(datExpr$X1, p = 0.75, list = FALSE)
train<- datExpr[inTraining,]
test<-datExpr[-inTraining,]
setLabels = c("Train", "Test");
multiExpr = list(Train = list(data = train), Test = list(data = test));
multiColor = list(Train =moduleColors  );
nSets = 2
mp = modulePreservation(multiExpr, multiColor,
                        referenceNetworks = 1,
                        nPermutations = 20,
                        randomSeed = 1,
                        quickCor = 0,
                        verbose = 3)
##save(mp,file = "mp.Rda") 这个很费时间，建议保存一下，nPermutations官网上给了200，我为节省时间给了20
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )


# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}






samples=read.csv('Sam_info.txt',sep = '\t',row.names = 1)

## 5 临床性状和表达热图

##表型与模块相关性##
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.9,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))













