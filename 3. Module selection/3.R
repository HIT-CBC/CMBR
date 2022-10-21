#模块选择
load("F:/002/code/data/methy_normal.RData")
load("F:/002/code/data/her2n/methy_her2n.RData")
load("F:/002/code/data/her2p/methy_her2p.RData")

test = cbind(methy_her2p,methy_normal)

Pvalue = apply(test,1,function(x){
  t.test(x[1:ncol(methy_her2p)],x[(ncol(methy_her2p)+1):(ncol(methy_her2p)+ncol(methy_normal))])$p.value
})
gap = apply(methy_her2p,1,mean)-apply(methy_normal,1,mean)

fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))

her2p_cg = rownames(methy_her2p)[which(fdr <0.05 & abs(gap) >=0.2)]
her2p_ttest = data.frame(Pvalue,fdr,gap)
write.csv(her2p_ttest,file = "F:/002/code/data/her2p_ttest.csv")

test = cbind(methy_her2n,methy_normal)

Pvalue = apply(test,1,function(x){
  t.test(x[1:ncol(methy_her2n)],x[(ncol(methy_her2n)+1):(ncol(methy_her2n)+ncol(methy_normal))])$p.value
})
gap = apply(methy_her2n,1,mean)-apply(methy_normal,1,mean)

fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))

her2n_cg = rownames(methy_her2n)[which(fdr <0.05 & abs(gap) >=0.2)]
her2n_ttest = data.frame(Pvalue,fdr,gap)
write.csv(her2n_ttest,file = "F:/002/code/data/her2n_ttest.csv")

hrp_cg = union(her2n_cg,her2p_cg)


#WGCNA
cgtest = read.table(file ="F:/002//002/data/cgtest.txt",header = F,sep = "\t")
load("F:/002/code/data/methy_cancer.RData")


cg_gene = cgtest$V1[!(grepl(";",cgtest$V3) | cgtest$V3 == "" )]
cg = intersect(cg_gene,hrp_cg)

hr_cg = methy_cancer[cg,]
hr_cg$cg = rownames(hr_cg)

library(dplyr)

hr_gene = merge(hr_cg,cgtest,by.x = "cg",by.y = "V1",all.x = F,all.y = F)
hr_gene = hr_gene[,-c(which(colnames(hr_gene) %in% c("V2","V4","V5")))]

hr_gene =  as.data.frame(hr_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

hr_gene = hr_gene[,-2]

rownames(hr_gene) = hr_gene[,1]
hr_gene = hr_gene[,-1]

write.table(hr_gene,file = "F:/002/code/data/hr_gene.txt",quote = F,sep = "\t")

library(WGCNA)
library(stringr)
library(reshape2)

options(stringsAsFactors = FALSE)

hr_gene = read.table(file = "F:/002/code/data/hr_gene.txt",header = T,sep = "\t")
datExpr <- t(hr_gene)
datExpr = as.data.frame(datExpr)
##软阈值筛选##
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
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

##一步法网络构建：One-step network construction and module detection##
net = blockwiseModules(datExpr, power = 9, maxBlockSize = 8000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)

##绘画结果展示##
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

##结果保存
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "F:/002/code/data/AS-green-FPKM-02-networkConstruction-auto.RData")
write(colnames(datExpr),file = "F:/002/code/data/gene_list.txt",sep = "\n")
######david,kegg通路分析

kegg = read.table(file = "F:/002/code/data/kegg.txt",sep = "\t",header = T,quote  = "")
kegg = kegg[kegg$PValue <= 0.05,]

kegg_gene = unique(unlist(strsplit(kegg[,"Genes"],', ')))

gene_list = cbind(moduleLabels,moduleColors)

blue = rownames(gene_list)[which(gene_list[,"moduleColors"] == "blue")]
brown = rownames(gene_list)[which(gene_list[,"moduleColors"] == "brown")]
grey = rownames(gene_list)[which(gene_list[,"moduleColors"] == "grey")]
turquoise = rownames(gene_list)[which(gene_list[,"moduleColors"] == "turquoise")]
yellow = rownames(gene_list)[which(gene_list[,"moduleColors"] == "yellow")]

a = length(intersect(blue,kegg_gene))/length(blue)
b = length(intersect(brown,kegg_gene))/length(brown)
c = length(intersect(grey,kegg_gene))/length(grey)
d = length(intersect(turquoise,kegg_gene))/length(turquoise)
e = length(intersect(yellow,kegg_gene))/length(yellow)

#选择棕色模块
write(brown,file = "F:/002/code/data/brown_gene.txt",sep = "\n")
#############################
