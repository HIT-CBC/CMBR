######################################################################
#分亚组
#her2p
library(survival)
library(survminer)
library(philentropy)
library(proxy)
load("F:/002/code/data/her2p/her2p_gene.RData")
load("F:/002/code/data/her2p/her2p_geoc_gene.RData")
clinic_BRCA = read.table(file = "F:/002/code/data/clin_BRCA.txt",sep = "\t")

clinic_BRCA$PFS = as.numeric(clinic_BRCA$PFS)
clinic_BRCA$PFS.time = as.numeric(clinic_BRCA$PFS.time)
colnames(her2p_gene) = substr(colnames(her2p_gene),1,12)
clinic_BRCA$bcr_patient_barcode = gsub("-",".",clinic_BRCA$bcr_patient_barcode)

mygene = read.table(file = "F:/002/code/data/her2p/mygene.txt")
mygene = mygene[,1]

#层次聚类
surviavl_data = t(her2p_gene[rownames(her2p_gene)%in%mygene,])
dis_bray <- dist(surviavl_data
                 ,method = 'chebyshev'
)
clust_average <- hclust(dis_bray, method = 'ward.D2')
clust_average_cut <- cutree(clust_average, k = 3)
write.table(clust_average_cut,file = "F:/002/code/data/her2p/cluster_average_cut.txt",quote = F,sep = "\t")
save(dis_bray,file = "F:/002/code/data/her2p/dis_bray.RData")

#生存曲线
surviavl_data = as.data.frame(surviavl_data)
surviavl_data=cbind(surviavl_data,clust_average_cut)
surviavl_data$bcr_patient_barcode = rownames(surviavl_data)

surviavl_data = merge(surviavl_data,clinic_BRCA,by = "bcr_patient_barcode",all.x = T,all.y = F)

surv1=Surv(surviavl_data$PFS.time,surviavl_data$PFS)

fit=survfit(surv1~clust_average_cut,data=surviavl_data) 
ggsurvplot( fit, 
            data = surviavl_data,
            size = 1, # change line size 
            palette = "aaas",
            #conf.int = TRUE, # Add confidence interval 
            pval = TRUE, # Add p-value 
            risk.table = TRUE, # Add risk table 
            risk.table.col = "clust_average_cut",# Risk table color by groups 
            legend.labs = c(paste0("cluster",1:3)), # Change legend labels
            risk.table.height = 0.25, # Useful to change when you have multiple groups 
            ggtheme = theme_bw() 
)

#geo聚类
surviavl_data = t(her2p_geoc_gene[rownames(her2p_geoc_gene)%in%mygene,])

dis_bray <- dist(surviavl_data
                 ,method = 'chebyshev'
)
clust_average <- hclust(dis_bray, method = 'ward.D2')
clust_average_cut <- cutree(clust_average, k = 3)
write.table(clust_average_cut,file = "F:/002/code/data/her2p/geoc_cluster_average_cut.txt",quote = F,sep = "\t")
save(dis_bray,file = "F:/002/code/data/her2p/geoc_dis_bray.RData")

#tsne
library(Rtsne)
#tcga
group=read.table("F:/002/code/data/her2p/cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2p/dis_bray.RData")
cluster = group$x

tsne_obj <- Rtsne(dis_bray, is_distance = TRUE,perplexity = 4)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(group$x),name = rownames(group))
ggplot(tsne_data,aes(x = X, y = Y,color = cluster)) +
  geom_point()+
  scale_color_manual(values = c( "#EE0000FF","#3B4992FF","#008B4599"))

#geo
group=read.table("F:/002/code/data/her2p/geoc_cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2p/geoc_dis_bray.RData")
cluster = group$x

tsne_obj <- Rtsne(dis_bray, is_distance = TRUE,perplexity = 4)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(group$x),name = rownames(group))
ggplot(tsne_data,aes(x = X, y = Y,color = cluster)) +
  geom_point()+
  scale_color_manual(values = c( "#EE0000FF","#3B4992FF","#008B4599"))

##her2n
library(survival)
library(survminer)
library(philentropy)
library(proxy)
load("F:/002/code/data/her2n/her2n_gene.RData")
load("F:/002/code/data/her2n/her2n_geoc_gene.RData")
clinic_BRCA = read.table(file = "F:/002/code/data/clin_BRCA.txt",sep = "\t")

clinic_BRCA$PFS = as.numeric(clinic_BRCA$PFS)
clinic_BRCA$PFS.time = as.numeric(clinic_BRCA$PFS.time)
colnames(her2n_gene) = substr(colnames(her2n_gene),1,12)
clinic_BRCA$bcr_patient_barcode = gsub("-",".",clinic_BRCA$bcr_patient_barcode)

mygene = read.table(file = "F:/002/code/data/her2n/mygene.txt")
mygene = mygene[,1]

#层次聚类
surviavl_data = t(her2n_gene[rownames(her2n_gene)%in%mygene,])
dis_bray <- dist(surviavl_data
                 ,method = 'chebyshev'
)
clust_average <- hclust(dis_bray, method = 'ward.D2')
clust_average_cut <- cutree(clust_average, k = 2)
write.table(clust_average_cut,file = "F:/002/code/data/her2n/cluster_average_cut.txt",quote = F,sep = "\t")
save(dis_bray,file = "F:/002/code/data/her2n/dis_bray.RData")

#生存曲线
surviavl_data = as.data.frame(surviavl_data)
surviavl_data=cbind(surviavl_data,clust_average_cut)
surviavl_data$bcr_patient_barcode = rownames(surviavl_data)

surviavl_data = merge(surviavl_data,clinic_BRCA,by = "bcr_patient_barcode",all.x = T,all.y = F)

surv1=Surv(surviavl_data$PFS.time,surviavl_data$PFS)

fit=survfit(surv1~clust_average_cut,data=surviavl_data) 
ggsurvplot( fit, 
            data = surviavl_data,
            size = 1, # change line size 
            palette = "aaas",
            #conf.int = TRUE, # Add confidence interval 
            pval = TRUE, # Add p-value 
            risk.table = TRUE, # Add risk table 
            risk.table.col = "clust_average_cut",# Risk table color by groups 
            legend.labs = c(paste0("cluster",1:2)), # Change legend labels
            risk.table.height = 0.25, # Useful to change when you have multiple groups 
            ggtheme = theme_bw() 
)

#geo聚类
surviavl_data = t(her2n_geoc_gene[rownames(her2n_geoc_gene)%in%mygene,])

dis_bray <- dist(surviavl_data
                 ,method = 'chebyshev'
)
clust_average <- hclust(dis_bray, method = 'ward.D2')
clust_average_cut <- cutree(clust_average, k = 2)
write.table(clust_average_cut,file = "F:/002/code/data/her2n/geoc_cluster_average_cut.txt",quote = F,sep = "\t")
save(dis_bray,file = "F:/002/code/data/her2n/geoc_dis_bray.RData")

#tsne
library(Rtsne)
#tcga
group=read.table("F:/002/code/data/her2n/cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2n/dis_bray.RData")
cluster = group$x

tsne_obj <- Rtsne(dis_bray, is_distance = TRUE,perplexity = 4)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(group$x),name = rownames(group))
ggplot(tsne_data,aes(x = X, y = Y,color = cluster)) +
  geom_point()+
  scale_color_manual(values = c( "#EE0000FF","#3B4992FF","#008B4599"))

#geo
group=read.table("F:/002/code/data/her2n/geoc_cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2n/geoc_dis_bray.RData")
cluster = group$x

tsne_obj <- Rtsne(dis_bray, is_distance = TRUE,perplexity = 4)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(group$x),name = rownames(group))
ggplot(tsne_data,aes(x = X, y = Y,color = cluster)) +
  geom_point()+
  scale_color_manual(values = c( "#EE0000FF","#3B4992FF","#008B4599"))

###########################################
#geo tcga 去批次做相关性
library(sva)
library(dplyr)

cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")

#her2p
load("F:/002/code/data/her2p/her2p_geoc_cg.RData")
load("F:/002/code/data/her2p/methy_her2p.RData")

cg = intersect(rownames(her2p_geoc_cg),rownames(methy_her2p))
methy_her2p = methy_her2p[cg,]
her2p_geoc_cg = her2p_geoc_cg[cg,]

methy = cbind(methy_her2p,her2p_geoc_cg)
sample = colnames(methy)
gse = c(rep("tcga",ncol(methy_her2p)),rep("geo",ncol(her2p_geoc_cg)))
combat_data =ComBat(dat=methy, batch=gse,par.prior=TRUE,prior.plots=FALSE)#校正批次效应
her2p_gene = as.data.frame(combat_data)

her2p_gene$cg = rownames(her2p_gene)
her2p_gene = merge(cgtest,her2p_gene,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2p_gene = her2p_gene[!(grepl(";",her2p_gene$V3) | her2p_gene$V3 == "" ),]
her2p_gene =  as.data.frame(her2p_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2p_gene = her2p_gene[,-(which(colnames(her2p_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2p_gene) = her2p_gene[,1]
her2p_gene = her2p_gene[,-1]

mygene = read.table(file = "F:/002/code/data/her2p/mygene.txt")
mygene = mygene[,1]
her2p_mygene = her2p_gene[rownames(her2p_gene) %in% mygene,]
her2p_mygene = as.data.frame(t(her2p_mygene))
her2p_mygene$patient = substr(rownames(her2p_mygene),1,12)

group1=read.table("F:/002/code/data/her2p/cluster_average_cut.txt",sep='\t',header=TRUE)
group2=read.table("F:/002/code/data/her2p/geoc_cluster_average_cut.txt",sep='\t',header=TRUE)
group1$patient = rownames(group1)
group2$patient = rownames(group2)
group = rbind(group1,group2)

test = merge(group,her2p_mygene,by = "patient")

mygene_tcga = test[grep("TCGA",test$patient),]
mygene_geo = test[grep("GSM",test$patient),]

ks.test(mygene_tcga$ACOT11,mygene_geo$ACOT11)

Pvalue = apply(test[,3:ncol(test)],2,function(x){
  t.test(x[1:nrow(mygene_geo)],x[(nrow(mygene_geo)+1):(nrow(mygene_geo)+nrow(mygene_tcga))])$p.value
})

#her2n
load("F:/002/code/data/her2n/her2n_geoc_cg.RData")
load("F:/002/code/data/her2n/methy_her2n.RData")

cg = intersect(rownames(her2n_geoc_cg),rownames(methy_her2n))
methy_her2n = methy_her2n[cg,]
her2n_geoc_cg = her2n_geoc_cg[cg,]

methy = cbind(methy_her2n,her2n_geoc_cg)
sample = colnames(methy)
gse = c(rep("tcga",ncol(methy_her2n)),rep("geo",ncol(her2n_geoc_cg)))
combat_data =ComBat(dat=methy, batch=gse,par.prior=TRUE,prior.plots=FALSE)#校正批次效应
her2n_gene = as.data.frame(combat_data)

her2n_gene$cg = rownames(her2n_gene)
her2n_gene = merge(cgtest,her2n_gene,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2n_gene = her2n_gene[!(grepl(";",her2n_gene$V3) | her2n_gene$V3 == "" ),]
her2n_gene =  as.data.frame(her2n_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2n_gene = her2n_gene[,-(which(colnames(her2n_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2n_gene) = her2n_gene[,1]
her2n_gene = her2n_gene[,-1]

mygene = read.table(file = "F:/002/code/data/her2n/mygene.txt")
mygene = mygene[,1]
her2n_mygene = her2n_gene[rownames(her2n_gene) %in% mygene,]
her2n_mygene = as.data.frame(t(her2n_mygene))
her2n_mygene$patient = substr(rownames(her2n_mygene),1,12)

group1=read.table("F:/002/code/data/her2n/cluster_average_cut.txt",sep='\t',header=TRUE)
group2=read.table("F:/002/code/data/her2n/geoc_cluster_average_cut.txt",sep='\t',header=TRUE)
group1$patient = rownames(group1)
group2$patient = rownames(group2)
group = rbind(group1,group2)

test = merge(group,her2n_mygene,by = "patient")

mygene_tcga = test[grep("TCGA",test$patient),]
mygene_geo = test[grep("GSM",test$patient),]

ks.test(mygene_tcga$ACOT11,mygene_geo$ACOT11)

Pvalue = apply(test[,3:ncol(test)],2,function(x){
  t.test(x[1:nrow(mygene_geo)],x[(nrow(mygene_geo)+1):(nrow(mygene_geo)+nrow(mygene_tcga))])$p.value
})

########################################
#免疫基因
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)

immgene = read.table(file = "F:/002/code/data/immGeneList.txt",sep = "\t",header = T,quote = "\"")

#her2p

group=read.table("F:/002/code/data/her2p/cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2p/her2p_gene.RData")

colnames(her2p_gene) = substr(colnames(her2p_gene),1,12)

test = her2p_gene[rownames(her2p_gene)%in%immgene$Symbol,]

test1 = test[,colnames(test)%in%(rownames(group)[group$x ==1])]
test2 = test[,colnames(test)%in%(rownames(group)[group$x ==2])]
test3 = test[,colnames(test)%in%(rownames(group)[group$x ==3])]

test$cluster1 = apply(test1,1,mean)
test$cluster2 = apply(test2,1,mean)
test$cluster3 = apply(test3,1,mean)

test$gene = rownames(test)

immdata = test[,c("gene","cluster1","cluster2","cluster3")]
write.table(immdata,file = "F:/002/code/data/her2p/immdata.txt",sep = "\t",quote = F,row.names = F)

#QDMR
immqdmr = read.table(file = "F:/002/code/data/her2p/QDMR/SpecificityTable.txt",sep = "\t",quote = "",header = T)
genelist = immqdmr$gene
immtable = her2p_gene[rownames(her2p_gene)%in%genelist,]
immtable = as.data.frame(t(immtable))
if(all(rownames(group) == rownames(immtable))){
  immtable$labels = as.character(group$x)
}


mydata<-melt(immtable,id.vars='labels')
colnames(mydata) = c("cluster","immgene","value")

mydata %>%
  ggplot(aes(x = immgene,y = value,fill = cluster))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "value")+
  scale_x_discrete(name = "immgene") +
  ggtitle("Comparison of immune cell infiltration") +coord_flip()+
  theme_bw()+scale_fill_manual(values =c( "#684e94","#6a855b","#5091c0"))+
  stat_compare_means(method = "anova",
                     #label = "p.format",##星号设置
                     label = "p.signif",
                     hide.ns=TRUE)+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold")) 

#her2n
group=read.table("F:/002/code/data/her2n/cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2n/her2n_gene.RData")

colnames(her2n_gene) = substr(colnames(her2n_gene),1,12)

test = her2n_gene[rownames(her2n_gene)%in%immgene$Symbol,]

test1 = test[,colnames(test)%in%(rownames(group)[group$x ==1])]
test2 = test[,colnames(test)%in%(rownames(group)[group$x ==2])]

test$cluster1 = apply(test1,1,mean)
test$cluster2 = apply(test2,1,mean)

test$gene = rownames(test)

immdata = test[,c("gene","cluster1","cluster2")]
write.table(immdata,file = "F:/002/code/data/her2n/immdata.txt",sep = "\t",quote = F,row.names = F)

#QDMR
immqdmr = read.table(file = "F:/002/code/data/her2n/QDMR/SpecificityTable.txt",sep = "\t",quote = "",header = T)
genelist = immqdmr$gene
immtable = her2n_gene[rownames(her2n_gene)%in%genelist,]
immtable = as.data.frame(t(immtable))
if(all(rownames(group) == rownames(immtable))){
  immtable$labels = as.character(group$x)
}


mydata<-melt(immtable,id.vars='labels')
colnames(mydata) = c("cluster","immgene","value")

mydata %>%
  ggplot(aes(x = immgene,y = value,fill = cluster))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "value")+
  scale_x_discrete(name = "immgene") +
  ggtitle("Comparison of immune cell infiltration") +coord_flip()+
  theme_bw()+scale_fill_manual(values =c( "#F39B7F","#B893CC"))+
  stat_compare_means(method = "anova",
                     #label = "p.format",##星号设置
                     label = "p.signif",
                     hide.ns=TRUE)+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold")) 
############################################
#fpkm2tpm
library(clusterProfiler)
library(org.Hs.eg.db)
library(DMwR2)

fpkm = read.table(file = "F:/002/code/data/TCGA-BRCA.htseq_fpkm.tsv",header=TRUE,sep='\t')
fpkm[,2:ncol(fpkm)] = 2^fpkm[,2:ncol(fpkm)]-1

gene=fpkm[,1]
gene=gsub("\\..*","",gene)

eg <- bitr(gene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")

fpkm$Ensembl_ID = gene
test = merge(x = eg,y = fpkm,by.x = "ENSEMBL",by.y = "Ensembl_ID",all.x = T,all.y = F)
test[,3:ncol(test)] = scale(test[,3:ncol(test)],center = F)

FPKM2TPM = function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#her2p
library(dplyr)

group=read.table("F:/002/code/data/her2p/cluster_average_cut.txt",sep='\t',header=TRUE)

patient = rownames(group)
samples = colnames(test)
samples = samples[as.numeric(substr(samples,14,15)) < 10]
samples = samples[substr(samples,1,12) %in% patient]

t = which(substr(samples,1,12) %in% substr(samples[which(duplicated(substr(samples,1,12)))],1,12))
n = samples[t][order(samples[t])][(1:(length(samples[t])/2))*2]
samples = samples[-c(which(samples %in% n))]

her2p_fpkm = test[,c(2,which(colnames(test)%in%samples))]
her2p_fpkm = her2p_fpkm[rowSums(her2p_fpkm[,2:ncol(her2p_fpkm)] == 0) <= (ncol(her2p_fpkm)-1)*0.7,]

her2p_tpm <- apply(her2p_fpkm[,2:ncol(her2p_fpkm)],2,FPKM2TPM)
her2p_tpm = as.data.frame(her2p_tpm)
her2p_tpm$SYMBOL = her2p_fpkm$SYMBOL

her2p_tpm =  her2p_tpm %>% group_by(SYMBOL) %>% summarise_each(funs(mean))

her2p_tpm = as.data.frame(her2p_tpm)
rownames(her2p_tpm) = her2p_tpm$SYMBOL
her2p_tpm = her2p_tpm[,-1]


save(her2p_tpm,file = "F:/002/code/data/her2p/her2p_tpm.RData")

#her2n
library(dplyr)

group=read.table(paste0("F:/002/code/data/her2n/cluster_average_cut.txt"),sep='\t',header=TRUE)

patient = rownames(group)
samples = colnames(test)
samples = samples[as.numeric(substr(samples,14,15)) < 10]
samples = samples[substr(samples,1,12) %in% patient]

t = which(substr(samples,1,12) %in% substr(samples[which(duplicated(substr(samples,1,12)))],1,12))
n = samples[t][order(samples[t])][(1:(length(samples[t])/2))*2]
samples = samples[-c(which(samples %in% n))]

her2n_fpkm = test[,c(2,which(colnames(test)%in%samples))]
her2n_fpkm = her2n_fpkm[rowSums(her2n_fpkm[,2:ncol(her2n_fpkm)] == 0) <= (ncol(her2n_fpkm)-1)*0.7,]

her2n_tpm <- apply(her2n_fpkm[,2:ncol(her2n_fpkm)],2,FPKM2TPM)
her2n_tpm = as.data.frame(her2n_tpm)
her2n_tpm$SYMBOL = her2n_fpkm$SYMBOL

her2n_tpm =  her2n_tpm %>% group_by(SYMBOL) %>% summarise_each(funs(mean))

her2n_tpm = as.data.frame(her2n_tpm)
rownames(her2n_tpm) = her2n_tpm$SYMBOL
her2n_tpm = her2n_tpm[,-1]

save(her2n_tpm,file = "F:/002/code/data/her2n/her2n_tpm.RData")

#CIBERSORT

#her2p

group=read.table(paste0("F:/002/code/data/her2p/cluster_average_cut.txt"),sep='\t',header=TRUE)
load("F:/002/code/data/her2p/her2p_tpm.RData")

tpm1 = her2p_tpm[,substr(colnames(her2p_tpm),1,12) %in% rownames(group)[group$x == 1]]
tpm2 = her2p_tpm[,substr(colnames(her2p_tpm),1,12) %in% rownames(group)[group$x == 2]]
tpm3 = her2p_tpm[,substr(colnames(her2p_tpm),1,12) %in% rownames(group)[group$x == 3]]

write.table(tpm1,file = "F:/002/code/data/her2p/tpm1.txt",sep = "\t")
write.table(tpm2,file = "F:/002/code/data/her2p/tpm2.txt",sep = "\t")
write.table(tpm3,file = "F:/002/code/data/her2p/tpm3.txt",sep = "\t")

setwd(dir = "F:/002/code/data/her2p/")

CIBERSORT("LM7.txt","tpm1.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_7.txt")
CIBERSORT("LM7.txt","tpm2.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_7.txt")
CIBERSORT("LM7.txt","tpm3.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_7.txt")


#her2n

group=read.table(paste0("F:/002/code/data/her2n/cluster_average_cut.txt"),sep='\t',header=TRUE)
load("F:/002/code/data/her2n/her2n_tpm.RData")

tpm1 = her2n_tpm[,substr(colnames(her2n_tpm),1,12) %in% rownames(group)[group$x == 1]]
tpm2 = her2n_tpm[,substr(colnames(her2n_tpm),1,12) %in% rownames(group)[group$x == 2]]

write.table(tpm1,file = "F:/002/code/data/her2n/tpm1.txt",sep = "\t")
write.table(tpm2,file = "F:/002/code/data/her2n/tpm2.txt",sep = "\t")

setwd(dir = "F:/002/code/data/her2n/")

CIBERSORT("LM7.txt","tpm1.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_7.txt")
CIBERSORT("LM7.txt","tpm2.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_7.txt")

###########################################
#geo  免疫基因
library(ggplot2)
library(ggpubr)

#her2p
group=read.table("F:/002/code/data/her2p/geoc_cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2p/her2p_geoc_gene.RData")

immqdmr = read.table(file = "F:/002/code/data/her2p/QDMR/SpecificityTable.txt",sep = "\t",quote = "",header = T)
genelist = immqdmr$gene
immtable = her2p_geoc_gene[rownames(her2p_geoc_gene)%in%genelist,]
immtable = as.data.frame(t(immtable))
if(all(rownames(group) == rownames(immtable))){
  immtable$labels = as.character(group$x)
}

e<-immtable %>% 
  dplyr::filter(labels %in% c("1","2","3")) %>% 
  ggviolin(x = "labels", y = c(colnames(immtable)[1:(ncol(immtable)-1)]),
           combine=T, select=c("1","2","3"),order=c("1","2","3"),
           palette =c( "#684e94","#6a855b","#5091c0"), color="labels",shape="labels",
           ylab="beta_value",xlab=FALSE,
           panel.labs = list(labels=c("1","2","3")),
           font.label = list(size = 14, face = "bold", color ="red"),
           add = "jitter", add.params = list(fill = "white"))
e+stat_compare_means(method = "kruskal.test",
                     #label = "p.format",##星号设置
                     label = "p.signif",
                     hide.ns = TRUE
)+theme_gray(base_size = 14)


#her2n
group=read.table("F:/002/code/data/her2n/geoc_cluster_average_cut.txt",sep='\t',header=TRUE)
load("F:/002/code/data/her2n/her2n_geoc_gene.RData")

immqdmr = read.table(file = "F:/002/code/data/her2n/QDMR/SpecificityTable.txt",sep = "\t",quote = "",header = T)
genelist = immqdmr$gene
immtable = her2n_geoc_gene[rownames(her2n_geoc_gene)%in%genelist,]
immtable = as.data.frame(t(immtable))
if(all(rownames(group) == rownames(immtable))){
  immtable$labels = as.character(group$x)
}

e<-immtable %>% 
  dplyr::filter(labels %in% c("1","2")) %>% 
  ggviolin(x = "labels", y = c(colnames(immtable)[1:(ncol(immtable)-1)]),
           combine=T, select=c("1","2"),order=c("1","2"),
           palette =c( "#684e94","#6a855b","#5091c0"), color="labels",shape="labels",
           ylab="beta_value",xlab=FALSE,
           panel.labs = list(labels=c("1","2")),
           font.label = list(size = 14, face = "bold", color ="red"),
           add = "jitter", add.params = list(fill = "white"))
e+stat_compare_means(method = "kruskal.test",
                     #label = "p.format",##星号设置
                     label = "p.signif",
                     hide.ns = TRUE
)+theme_gray(base_size = 14)

############################################
#免疫检查点
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)

#her2p
load("F:/002/code/data/her2p/her2p_tpm.RData")

group=read.table(paste0("F:/002/code/data/her2p/cluster_average_cut.txt"),sep='\t',header=TRUE)
group$sample_id = rownames(group)

colnames(her2p_tpm) = substr(colnames(her2p_tpm),1,12)

checkpoint = c("CTLA4","PDCD1","CD274","CD226","TIGIT","CRTAM",
               "CD96","CD200R1","KLRD1","KLRC1","TNFRSF9","LAYN","LILRB1")

checkpoint %in% rownames(her2p_tpm)

test = her2p_tpm[rownames(her2p_tpm) %in% checkpoint,]
test = as.data.frame(t(test))

test$sample_id = rownames(test)

test = merge(test,group,by = "sample_id",all.x = F,all.y = F)
rownames(test) = test$sample_id
test = test[,-1]
test$x = as.character(test$x)

her2p_checkpoint =  as.data.frame(test %>% group_by(x) %>% summarise_each(funs(median)))
her2p_checkpoint = her2p_checkpoint[,-1]
her2p_checkpoint = as.data.frame(t(her2p_checkpoint))
colnames(her2p_checkpoint) = paste0(rep("cluster",ncol(her2p_checkpoint)),1:ncol(her2p_checkpoint))
write.csv(her2p_checkpoint,file = "F:/002/code/data/her2p_checkpoint.csv")

mydata<-melt(test,id.vars='x')
colnames(mydata) = c("cluster","checkpoint","value")

mydata %>%
  ggplot(aes(x = checkpoint,y = value,fill = cluster))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "value")+
  scale_x_discrete(name = "checkpoint") +
  ggtitle("Comparison of immune cell infiltration") +coord_flip(ylim = c(0,30))+
  theme_bw()+scale_fill_manual(values =c( "#684e94","#6a855b","#5091c0"))+
  stat_compare_means(method = "kruskal.test",
                     #label = "p.format",##星号设置
                     label = "p.signif",label.y = 30,
                     hide.ns=TRUE       )+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold")) 

#her2n
load("F:/002/code/data/her2n/her2n_tpm.RData")

group=read.table(paste0("F:/002/code/data/her2n/cluster_average_cut.txt"),sep='\t',header=TRUE)
group$sample_id = rownames(group)

colnames(her2n_tpm) = substr(colnames(her2n_tpm),1,12)

checkpoint = c("CTLA4","PDCD1","CD274","CD226","TIGIT","CRTAM",
               "CD96","CD200R1","KLRD1","KLRC1","TNFRSF9","LAYN","LILRB1")

checkpoint %in% rownames(her2n_tpm)

test = her2n_tpm[rownames(her2n_tpm) %in% checkpoint,]
test = as.data.frame(t(test))

test$sample_id = rownames(test)

test = merge(test,group,by = "sample_id",all.x = F,all.y = F)
rownames(test) = test$sample_id
test = test[,-1]
test$x = as.character(test$x)

her2n_checkpoint =  as.data.frame(test %>% group_by(x) %>% summarise_each(funs(median)))
her2n_checkpoint = her2n_checkpoint[,-1]
her2n_checkpoint = as.data.frame(t(her2n_checkpoint))
colnames(her2n_checkpoint) = paste0(rep("cluster",ncol(her2n_checkpoint)),1:ncol(her2n_checkpoint))
write.csv(her2n_checkpoint,file = "F:/002/code/data/her2n_checkpoint.csv")

mydata<-melt(test,id.vars='x')
colnames(mydata) = c("cluster","checkpoint","value")

mydata %>%
  ggplot(aes(x = checkpoint,y = value,fill = cluster))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "value")+
  scale_x_discrete(name = "checkpoint") +
  ggtitle("Comparison of immune cell infiltration") +coord_flip(ylim = c(0,30))+
  theme_bw()+scale_fill_manual(values =c( "#F39B7F","#B893CC"))+
  stat_compare_means(method = "kruskal.test",
                     #label = "p.format",##星号设置
                     label = "p.signif",label.y = 30,
                     hide.ns=TRUE       )+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold")) 
