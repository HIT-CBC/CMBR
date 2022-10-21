#根据IHC划分样本

TCGA_cell = read.table(file = "F:/002/code/data/TCGA_cell_2015.tsv",sep = "\t",quote = "",header = T)
TCGA_nature = read.table(file = "F:/002/code/data/TCGA_nature_2012_BRCA.tsv",sep = "\t",quote = "",header = T)
TCGA_firehose = read.table(file = "F:/002/code/data/TCGA_Firehose_Legacy.tsv",sep = "\t",quote = "",header = T)

TCGA_cell = TCGA_cell[-c(which(substr(TCGA_cell$Sample.ID,14,16)!="01")),]
TCGA_firehose = TCGA_firehose[-c(which(substr(TCGA_firehose$Sample.ID,14,16)!="01")),]

TCGA_cell = TCGA_cell[,c("Patient.ID","ER.Status.By.IHC","PR.status.by.ihc","IHC.HER2")]
TCGA_nature = TCGA_nature[,c("Patient.ID","ER.Status","PR.Status","HER2.Status")]
TCGA_firehose = TCGA_firehose[,c("Patient.ID","ER.Status.By.IHC","PR.status.by.ihc","IHC.HER2")]

clin = read.table("F:/002/code/data/clin.txt",sep = "\t",quote = "",header = T)

subtype_ihc = clin[,c("bcr_patient_barcode","er_status_by_ihc","pr_status_by_ihc","her2_status_by_ihc")]

colnames(subtype_ihc) = c("ID","ER","PR","HER2")
colnames(TCGA_cell) = c("ID","ER","PR","HER2")
colnames(TCGA_nature) = c("ID","ER","PR","HER2")
colnames(TCGA_firehose) = c("ID","ER","PR","HER2")

a = merge(subtype_ihc,TCGA_cell,by = "ID",all.x = T,all.y = T)
b = merge(TCGA_nature,TCGA_firehose,by = "ID",all.x = T,all.y = T)

IHC_subtype = merge(a,b,by = "ID",all.x = T,all.y = T)
IHC_subtype = IHC_subtype[,c("ID","ER.x.x","ER.x.y","ER.y.x","ER.y.y","PR.x.x","PR.x.y","PR.y.x","PR.y.y","HER2.x.x","HER2.x.y","HER2.y.x","HER2.y.y")]


ER_P = (IHC_subtype$ER.x.x == "Positive" )|(IHC_subtype$ER.x.y == "Positive" )|(IHC_subtype$ER.y.x == "Positive" )|(IHC_subtype$ER.y.y == "Positive" )
ER_N = (IHC_subtype$ER.x.x == "Negative" )|(IHC_subtype$ER.x.y == "Negative" )|(IHC_subtype$ER.y.x == "Negative" )|(IHC_subtype$ER.y.y == "Negative" )

PR_P = (IHC_subtype$PR.x.x == "Positive" )|(IHC_subtype$PR.x.y == "Positive" )|(IHC_subtype$PR.y.x == "Positive" )|(IHC_subtype$PR.y.y == "Positive" )
PR_N = (IHC_subtype$PR.x.x == "Negative" )|(IHC_subtype$PR.x.y == "Negative" )|(IHC_subtype$PR.y.x == "Negative" )|(IHC_subtype$PR.y.y == "Negative" )

HER2_P = (IHC_subtype$HER2.x.x == "Positive" )|(IHC_subtype$HER2.x.y == "Positive" )|(IHC_subtype$HER2.y.x == "Positive" )|(IHC_subtype$HER2.y.y == "Positive" )
HER2_N = (IHC_subtype$HER2.x.x == "Negative" )|(IHC_subtype$HER2.x.y == "Negative" )|(IHC_subtype$HER2.y.x == "Negative" )|(IHC_subtype$HER2.y.y == "Negative" )


IHC_subtype$ER[ER_P] = "Positive"
IHC_subtype$ER[ER_N] = "Negative"

IHC_subtype$PR[PR_P] = "Positive"
IHC_subtype$PR[PR_N] = "Negative"

IHC_subtype$HER2[HER2_P] = "Positive"
IHC_subtype$HER2[HER2_N] = "Negative"

IHC_subtype = IHC_subtype[,c("ID","ER","PR","HER2")]

IHC_subtype$subtype[(IHC_subtype$ER=='Positive' | IHC_subtype$PR=='Positive') & IHC_subtype$HER2=='Negative'] = 'her2n'
IHC_subtype$subtype[(IHC_subtype$ER=='Positive' | IHC_subtype$PR=='Positive') & IHC_subtype$HER2=='Positive'] = 'her2p'
IHC_subtype$subtype[(IHC_subtype$ER=='Negative' & IHC_subtype$PR=='Negative') & IHC_subtype$HER2=='Positive'] = 'HER2'
IHC_subtype$subtype[(IHC_subtype$ER=='Negative' & IHC_subtype$PR=='Negative') & IHC_subtype$HER2=='Negative'] = 'TNBC'

IHC_subtype$subtype[is.na(IHC_subtype$subtype)] = "unknow"
IHC_subtype$ID = gsub("-",".",IHC_subtype$ID)

save(IHC_subtype,file = "F:/002/code/data/IHC_subtype.RData")

###########################
#27k和450k去批次
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")

methy450 = read.table(file = "F:/002/code/data/TCGA-BRCA.methylation450.tsv",header = T,sep = "\t",nrows = 10000,row.names = 1)
methy450 = methy450[rownames(methy450) %in% cgtest$V1,]
methy27 = read.table(file = "F:/002/code/data/TCGA-BRCA.methylation27.tsv",header = T,sep = "\t",row.names = 1)

save(methy27,file = "F:/002/code/data/methy27.RData")
save(methy450,file = "F:/002/code/data/methy450.RData")

methy27 = methy27[rowSums(is.na(methy27)) <= (ncol(methy27)*0.7),]
methy450 = methy450[rowSums(is.na(methy450)) <= (ncol(methy450)*0.7),]

cg = intersect(rownames(methy27),rownames(methy450))
methy27 = methy27[cg,]
methy450 = methy450[cg,]

methy27_cancer = methy27[,as.numeric(substr(colnames(methy27),14,15)) < 10 ]
methy27_normal = methy27[,as.numeric(substr(colnames(methy27),14,15)) >= 10 ]

methy450_cancer = methy450[,as.numeric(substr(colnames(methy450),14,15)) < 10 ]
methy450_normal = methy450[,as.numeric(substr(colnames(methy450),14,15)) >= 10 ]


methy27_cancer_patient_id = substr(colnames(methy27_cancer),1,12)
methy450_cancer_patient_id = substr(colnames(methy450_cancer),1,12)

a=0
for (i in methy27_cancer_patient_id[duplicated(methy27_cancer_patient_id)]) {
  a = c(a,grep(i,colnames(methy27_cancer)))
}
a = a[-1]
colnames(methy27_cancer)[a]
methy27_cancer = methy27_cancer[,-a[2:3]]
a = 0
for (i in methy450_cancer_patient_id[duplicated(methy450_cancer_patient_id)]) {
  a = c(a,grep(i,colnames(methy450_cancer)))
}
a = a[-1]
colnames(methy450_cancer)[a]
methy450_cancer = methy450_cancer[,-a[c(1,3,6,8,10,12,14,16,18,19,22)]]

RNAseq_id = read.table(file = "F:/002/code/data/RNAseq_id.txt")
RNAseq_id$id = gsub("-",".",substr(RNAseq_id$V1,1,12))

methy450_cancer = methy450_cancer[,substr(colnames(methy450_cancer),1,12) %in% RNAseq_id$id]
methy27_cancer = methy27_cancer[,substr(colnames(methy27_cancer),1,12) %in% RNAseq_id$id]

library(DMwR2)
methy450_cancer <- knnImputation(methy450_cancer, k=10, scale = T, meth = "weighAvg")
methy27_cancer <- knnImputation(methy27_cancer, k=10, scale = T, meth = "weighAvg")
methy450_normal <- knnImputation(methy450_normal, k=10, scale = T, meth = "weighAvg")
methy27_normal <- knnImputation(methy27_normal, k=10, scale = T, meth = "weighAvg")


methy = cbind(methy450_cancer,methy27_cancer,methy450_normal,methy27_normal)
sample = colnames(methy)
gse = c(rep("450k",ncol(methy450_cancer)),rep("27k",ncol(methy27_cancer)),rep("450k",ncol(methy450_normal)),rep("27k",ncol(methy27_normal)))
class = c(rep("cancer",ncol(methy450_cancer)),rep("cancer",ncol(methy27_cancer)),rep("normal",ncol(methy450_normal)),rep("normal",ncol(methy27_normal)))
bdata = as.data.frame(cbind(sample,gse,class))
mod = model.matrix(~as.factor(class), data=bdata)#根据表型文件建立模型
library(sva)
combat_data =ComBat(dat=methy, batch=bdata$gse, mod=mod,par.prior=TRUE,prior.plots=FALSE)#校正批次效应

save(combat_data,file = "F:/002/code/data/combat_data.RData")

#####################################################################
#划分样本
load("F:/002/code/data/IHC_subtype.RData")
load("F:/002/code/data/combat_data.RData")

normal_id = colnames(combat_data)[as.numeric(substr(colnames(combat_data),14,15))>=10]
cancer_id = colnames(combat_data)[as.numeric(substr(colnames(combat_data),14,15))<10]
her2n_id = IHC_subtype$ID[IHC_subtype$subtype == "her2n"]
her2p_id = IHC_subtype$ID[IHC_subtype$subtype == "her2p"]
other_id = IHC_subtype$ID[IHC_subtype$subtype == "HER2" | IHC_subtype$subtype == "TNBC"]
her2_id = IHC_subtype$ID[IHC_subtype$subtype == "HER2" ]
tnbc_id = IHC_subtype$ID[ IHC_subtype$subtype == "TNBC"]

cancer_id = cancer_id[substr(cancer_id,1,12)%in%c(her2n_id,her2p_id,other_id)]
hrn_id = cancer_id[substr(cancer_id,1,12)%in%c(other_id)]
hrp_id = cancer_id[substr(cancer_id,1,12)%in%c(her2n_id,her2p_id)]
her2p_id= cancer_id[substr(cancer_id,1,12)%in%c(her2p_id)]
her2n_id= cancer_id[substr(cancer_id,1,12)%in%c(her2n_id)]

methy_cancer = as.data.frame(combat_data[,cancer_id])
methy_normal = as.data.frame(combat_data[,normal_id])
methy_hrp = as.data.frame(combat_data[,hrp_id])
methy_hrn = as.data.frame(combat_data[,hrn_id])
methy_her2p = as.data.frame(combat_data[,her2p_id])
methy_her2n = as.data.frame(combat_data[,her2n_id])


save(methy_cancer,file = "F:/002/code/data/methy_cancer.RData")
save(methy_normal,file = "F:/002/code/data/methy_normal.RData")
save(methy_hrp,file = "F:/002/code/data/methy_hrp.RData")
save(methy_hrn,file = "F:/002/code/data/methy_hrn.RData")

save(methy_her2p,file = "F:/002/code/data/her2p/methy_her2p.RData")
save(methy_her2n,file = "F:/002/code/data/her2n/methy_her2n.RData")

cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")

library(dplyr)
methy_her2p$cg = rownames(methy_her2p)
her2p_gene = merge(cgtest,methy_her2p,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2p_gene = her2p_gene[!(grepl(";",her2p_gene$V3) | her2p_gene$V3 == "" ),]
her2p_gene =  as.data.frame(her2p_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2p_gene = her2p_gene[,-(which(colnames(her2p_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2p_gene) = her2p_gene[,1]
her2p_gene = her2p_gene[,-1]
save(her2p_gene,file = "F:/002/code/data/her2p/her2p_gene.RData")


methy_her2n$cg = rownames(methy_her2n)
her2n_gene = merge(cgtest,methy_her2n,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2n_gene = her2n_gene[!(grepl(";",her2n_gene$V3) | her2n_gene$V3 == "" ),]
her2n_gene =  as.data.frame(her2n_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2n_gene = her2n_gene[,-(which(colnames(her2n_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2n_gene) = her2n_gene[,1]
her2n_gene = her2n_gene[,-1]
save(her2n_gene,file = "F:/002/code/data/her2n/her2n_gene.RData")

###############################################
#geo数据整合
#geo数据集1划分
methy_geo = read.table(file = "F:/002/code/data/geo_methy.txt",header = T,row.names = 1,sep = "\t")
clin_geo = read.table(file = "F:/002/code/data/geo_clin.txt",sep = "\t")
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")
load("F:/002/code/data/her2p/methy_her2p.RData")

library(dplyr)
library(DMwR2)

methy_geo = methy_geo[rownames(methy_geo)%in%rownames(methy_her2p),]

clin_geo = clin_geo[,-1]


clin_Cohort1 = clin_geo[,1:118]
clin_Cohort2 = clin_geo[,119:237]

clin_Cohort1 = as.data.frame(t(clin_Cohort1))
clin_Cohort2 = as.data.frame(t(clin_Cohort2))


clin_Cohort1 = clin_Cohort1[,c("V2","V14","V15","V22","V23")]
clin_Cohort2 = clin_Cohort2[,c("V2","V13","V24","V27","V28")]


clin_Cohort2$subtype[(clin_Cohort2$V24=='er_ihc_status: 1' | clin_Cohort2$V27=='pgr_status: 1') & clin_Cohort2$V28=='her2_status: 0'] = 'her2n'
clin_Cohort2$subtype[(clin_Cohort2$V24=='er_ihc_status: 1' | clin_Cohort2$V27=='pgr_status: 1') & clin_Cohort2$V28=='her2_status: 1'] = 'her2p'
clin_Cohort2$subtype[(clin_Cohort2$V24=='er_ihc_status: 0' & clin_Cohort2$V27=='pgr_status: 0') & clin_Cohort2$V28=='her2_status: 1'] = 'HER2'
clin_Cohort2$subtype[(clin_Cohort2$V24=='er_ihc_status: 0' & clin_Cohort2$V27=='pgr_status: 0') & clin_Cohort2$V28=='her2_status: 0'] = 'TNBC'


clin_Cohort1$subtype[clin_Cohort1$V22=='er_ihc_status: 1' & clin_Cohort1$V23 =='her2_status: 0'] = 'her2n'
clin_Cohort1$subtype[clin_Cohort1$V22=='er_ihc_status: 1' & clin_Cohort1$V23 =='her2_status: 1'] = 'her2p'


geo_her2n = c(clin_Cohort1$V2[clin_Cohort1$subtype == "her2n"],clin_Cohort2$V2[clin_Cohort2$subtype == "her2n"])
geo_her2p = c(clin_Cohort1$V2[clin_Cohort1$subtype == "her2p"],clin_Cohort2$V2[clin_Cohort2$subtype == "her2p"])

geo_her2n = na.omit(geo_her2n)
geo_her2p = na.omit(geo_her2p)

methy_geo_her2n = methy_geo[,geo_her2n]
methy_geo_her2p = methy_geo[,geo_her2p]

methy_geo_her2n = apply(methy_geo_her2n,2,function(x) {
  gsub("null",NA,x)
})
methy_geo_her2p = apply(methy_geo_her2p,2,function(x) {
  gsub("null",NA,x)
})

methy_geo_her2n = methy_geo_her2n[rowSums(is.na(methy_geo_her2n)) <= (ncol(methy_geo_her2n)*0.7),]
methy_geo_her2p = methy_geo_her2p[rowSums(is.na(methy_geo_her2p)) <= (ncol(methy_geo_her2p)*0.7),]

cg = rownames(methy_geo_her2n)

methy_geo_her2n = apply(methy_geo_her2n,2,as.numeric)
methy_geo_her2p = apply(methy_geo_her2p,2,as.numeric)

methy_geo_her2n = as.data.frame(methy_geo_her2n)
methy_geo_her2p = as.data.frame(methy_geo_her2p)

rownames(methy_geo_her2n) = cg
rownames(methy_geo_her2p) = cg

methy_geo_her2n <- knnImputation(methy_geo_her2n, k=10, scale = T, meth = "weighAvg")
methy_geo_her2p <- knnImputation(methy_geo_her2p, k=10, scale = T, meth = "weighAvg")
her2p_geo_cg = methy_geo_her2p
her2n_geo_cg = methy_geo_her2n
save(her2p_geo_cg,file = "F:/002/code/data/her2p/her2p_geo_cg.RData")
save(her2n_geo_cg,file = "F:/002/code/data/her2n/her2n_geo_cg.RData")
methy_geo_her2n$cg = rownames(methy_geo_her2n)
methy_geo_her2p$cg = rownames(methy_geo_her2p)

her2n_geo_gene = merge(cgtest,methy_geo_her2n,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2n_geo_gene = her2n_geo_gene[!(grepl(";",her2n_geo_gene$V3) | her2n_geo_gene$V3 == "" ),]
her2n_geo_gene =  as.data.frame(her2n_geo_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2n_geo_gene = her2n_geo_gene[,-(which(colnames(her2n_geo_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2n_geo_gene) = her2n_geo_gene[,1]
her2n_geo_gene = her2n_geo_gene[,-1]


her2p_geo_gene = merge(cgtest,methy_geo_her2p,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2p_geo_gene = her2p_geo_gene[!(grepl(";",her2p_geo_gene$V3) | her2p_geo_gene$V3 == "" ),]
her2p_geo_gene =  as.data.frame(her2p_geo_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2p_geo_gene = her2p_geo_gene[,-(which(colnames(her2p_geo_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2p_geo_gene) = her2p_geo_gene[,1]
her2p_geo_gene = her2p_geo_gene[,-1]


save(her2n_geo_gene,file = "F:/002/code/data/her2n/her2n_geo_gene.RData")
save(her2p_geo_gene,file = "F:/002/code/data/her2p/her2p_geo_gene.RData")

#geo2
methy_geo = read.table(file = "F:/002/code/data/geo2_methy.txt",header = T,row.names = 1,sep = "\t")
clin_geo = read.table(file = "F:/002/code/data/geo2_clin.txt",sep = "\t")
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")
load("F:/002/code/data/her2p/methy_her2p.RData")

library(dplyr)
library(DMwR2)

methy_geo = methy_geo[rownames(methy_geo)%in%rownames(methy_her2p),]

clin_geo = clin_geo[,-1]

clin_geo = as.data.frame(t(clin_geo))

clin_geo = clin_geo[,c("V2","V17","V18")]

clin_geo$subtype[clin_geo$V17=='er: 1' & clin_geo$V18 =='her2: 0'] = 'her2n'
clin_geo$subtype[clin_geo$V17=='er: 1' & clin_geo$V18 =='her2: 1'] = 'her2p'

her2n_geo2_gene = methy_geo[,colnames(methy_geo) %in% clin_geo$V2[clin_geo$subtype == 'her2n']]
her2n_geo2_cg = her2n_geo2_gene
save(her2n_geo2_cg,file = "F:/002/code/data/her2n/her2n_geo2_cg.RData")
her2n_geo2_gene$cg = rownames(her2n_geo2_gene)
her2n_geo2_gene = merge(cgtest,her2n_geo2_gene,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2n_geo2_gene = her2n_geo2_gene[!(grepl(";",her2n_geo2_gene$V3) | her2n_geo2_gene$V3 == "" ),]
her2n_geo2_gene =  as.data.frame(her2n_geo2_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2n_geo2_gene = her2n_geo2_gene[,-(which(colnames(her2n_geo2_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2n_geo2_gene) = her2n_geo2_gene[,1]
her2n_geo2_gene = her2n_geo2_gene[,-1]


her2p_geo2_gene = methy_geo[,colnames(methy_geo) %in% clin_geo$V2[clin_geo$subtype == 'her2p']]
her2p_geo2_cg = her2p_geo2_gene
save(her2p_geo2_cg,file = "F:/002/code/data/her2p/her2p_geo2_cg.RData")
her2p_geo2_gene$cg = rownames(her2p_geo2_gene)
her2p_geo2_gene = merge(cgtest,her2p_geo2_gene,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2p_geo2_gene = her2p_geo2_gene[!(grepl(";",her2p_geo2_gene$V3) | her2p_geo2_gene$V3 == "" ),]
her2p_geo2_gene =  as.data.frame(her2p_geo2_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2p_geo2_gene = her2p_geo2_gene[,-(which(colnames(her2p_geo2_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2p_geo2_gene) = her2p_geo2_gene[,1]
her2p_geo2_gene = her2p_geo2_gene[,-1]

save(her2n_geo2_gene,file = "F:/002/code/data/her2n/her2n_geo2_gene.RData")
save(her2p_geo2_gene,file = "F:/002/code/data/her2p/her2p_geo2_gene.RData")



#geo3
methy_geo = read.table(file = "F:/002/code/data/geo3_methy.txt",header = T,row.names = 1,sep = "\t")
clin_geo = read.table(file = "F:/002/code/data/geo3_clin.txt",sep = "\t")
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")
load("F:/002/code/data/her2p/methy_her2p.RData")

library(dplyr)
library(DMwR2)

methy_geo = methy_geo[rownames(methy_geo)%in%rownames(methy_her2p),]
t = rownames(methy_geo)
methy_geo = apply(methy_geo,2,as.numeric)
methy_geo = as.data.frame(methy_geo)
rownames(methy_geo) = t

clin_geo = clin_geo[,-1]

clin_geo = as.data.frame(t(clin_geo))

clin_geo = clin_geo[,c("V2","V12","V13","V14","V15","V16")]

clin_geo$subtype[(clin_geo$V12=='er: +' | clin_geo$V13=='pr: +') & clin_geo$V14=='her2: 0'] = 'her2n'
clin_geo$subtype[(clin_geo$V12=='er: +' | clin_geo$V13=='pr: +') & clin_geo$V14=='her2: +'] = 'her2p'

her2n_geo3_gene = methy_geo[,colnames(methy_geo) %in% clin_geo$V2[clin_geo$subtype == 'her2n']]
her2n_geo3_gene = knnImputation(her2n_geo3_gene, k=10, scale = F, meth = "weighAvg")
her2n_geo3_cg = her2n_geo3_gene
save(her2n_geo3_cg,file = "F:/002/code/data/her2n/her2n_geo3_cg.RData")
her2n_geo3_gene$cg = rownames(her2n_geo3_gene)
her2n_geo3_gene = merge(cgtest,her2n_geo3_gene,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2n_geo3_gene = her2n_geo3_gene[!(grepl(";",her2n_geo3_gene$V3) | her2n_geo3_gene$V3 == "" ),]
her2n_geo3_gene =  as.data.frame(her2n_geo3_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2n_geo3_gene = her2n_geo3_gene[,-(which(colnames(her2n_geo3_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2n_geo3_gene) = her2n_geo3_gene[,1]
her2n_geo3_gene = her2n_geo3_gene[,-1]


her2p_geo3_gene = methy_geo[,colnames(methy_geo) %in% clin_geo$V2[clin_geo$subtype == 'her2p']]
her2p_geo3_gene = knnImputation(her2p_geo3_gene, k=10, scale = F, meth = "weighAvg")
her2p_geo3_cg = her2p_geo3_gene
save(her2p_geo3_cg,file = "F:/002/code/data/her2p/her2p_geo3_cg.RData")
her2p_geo3_gene$cg = rownames(her2p_geo3_gene)
her2p_geo3_gene = merge(cgtest,her2p_geo3_gene,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2p_geo3_gene = her2p_geo3_gene[!(grepl(";",her2p_geo3_gene$V3) | her2p_geo3_gene$V3 == "" ),]
her2p_geo3_gene =  as.data.frame(her2p_geo3_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2p_geo3_gene = her2p_geo3_gene[,-(which(colnames(her2p_geo3_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2p_geo3_gene) = her2p_geo3_gene[,1]
her2p_geo3_gene = her2p_geo3_gene[,-1]

save(her2n_geo3_gene,file = "F:/002/code/data/her2n/her2n_geo3_gene.RData")
save(her2p_geo3_gene,file = "F:/002/code/data/her2p/her2p_geo3_gene.RData")


#geo4
methy_geo = read.table(file = "F:/002/code/data/geo4_methy.txt",header = T,row.names = 1,sep = "\t")
clin_geo = read.table(file = "F:/002/code/data/geo4_clin.txt",sep = "\t")
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")
load("F:/002/code/data/her2p/methy_her2p.RData")

library(dplyr)
library(DMwR2)

methy_geo = methy_geo[rownames(methy_geo)%in%rownames(methy_her2p),]
t = rownames(methy_geo)
methy_geo = apply(methy_geo,2,as.numeric)
methy_geo = as.data.frame(methy_geo)
rownames(methy_geo) = t

clin_geo = clin_geo[,-1]

clin_geo = as.data.frame(t(clin_geo))

clin_geo = clin_geo[,c("V2","V16","V17")]

clin_geo$subtype[(clin_geo$V16=='er: Positive' | clin_geo$V17=='pgr: Positive')] = 'her2p'

her2p_geo4_gene = methy_geo[,colnames(methy_geo) %in% clin_geo$V2[clin_geo$subtype == 'her2p']]
her2p_geo4_gene = knnImputation(her2p_geo4_gene, k=10, scale = F, meth = "weighAvg")
her2p_geo4_cg = her2p_geo4_gene
save(her2p_geo4_cg,file = "F:/002/code/data/her2p/her2p_geo4_cg.RData")
her2p_geo4_gene$cg = rownames(her2p_geo4_gene)
her2p_geo4_gene = merge(cgtest,her2p_geo4_gene,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2p_geo4_gene = her2p_geo4_gene[!(grepl(";",her2p_geo4_gene$V3) | her2p_geo4_gene$V3 == "" ),]
her2p_geo4_gene =  as.data.frame(her2p_geo4_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2p_geo4_gene = her2p_geo4_gene[,-(which(colnames(her2p_geo4_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2p_geo4_gene) = her2p_geo4_gene[,1]
her2p_geo4_gene = her2p_geo4_gene[,-1]

save(her2p_geo4_gene,file = "F:/002/code/data/her2p/her2p_geo4_gene.RData")


#合并geo数据
library(pheatmap)

#热图看差异程度
#her2p

load("F:/002/code/data/her2p/her2p_geo_cg.RData")
load("F:/002/code/data/her2p/her2p_geo2_cg.RData")
load("F:/002/code/data/her2p/her2p_geo3_cg.RData")
load("F:/002/code/data/her2p/her2p_geo4_cg.RData")

cg = intersect(intersect(rownames(her2p_geo_cg),rownames(her2p_geo2_cg)),
                 intersect(rownames(her2p_geo3_cg),rownames(her2p_geo4_cg)))

her2p_geo_cg = her2p_geo_cg[cg,]
her2p_geo2_cg = her2p_geo2_cg[cg,]
her2p_geo3_cg = her2p_geo3_cg[cg,]
her2p_geo4_cg = her2p_geo4_cg[cg,]

diff = cbind(her2p_geo_cg,her2p_geo2_cg,her2p_geo3_cg,her2p_geo4_cg)

anno_col=data.frame(sampleType=factor(rep(c("geo1","geo2","geo3","geo4"),
                                          c(ncol(her2p_geo_cg),ncol(her2p_geo2_cg),
                                            ncol(her2p_geo3_cg),ncol(her2p_geo4_cg)))))
rownames(anno_col)=colnames(diff)

pheatmap(diff[sample(1:nrow(diff),1000),],scale="none",annotation_col=anno_col,#annotation_colors=ann_color,
         clustering_method = 'ward.D2',
         show_colnames=F,show_rownames=F,
         fontsize_row=6,fontsize=10,
         legend_breaks=c(0.2,0.8),legend_labels=c('low','high'),
         cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c('#436eee','white','#EE0000'))(100)
)


#her2n

load("F:/002/code/data/her2n/her2n_geo_cg.RData")
load("F:/002/code/data/her2n/her2n_geo2_cg.RData")
load("F:/002/code/data/her2n/her2n_geo3_cg.RData")

cg = intersect(intersect(rownames(her2n_geo_cg),rownames(her2n_geo2_cg)),
                 rownames(her2n_geo3_cg))

her2n_geo_cg = her2n_geo_cg[cg,]
her2n_geo2_cg = her2n_geo2_cg[cg,]
her2n_geo3_cg = her2n_geo3_cg[cg,]

diff = cbind(her2n_geo_cg,her2n_geo2_cg,her2n_geo3_cg)

anno_col=data.frame(sampleType=factor(rep(c("geo1","geo2","geo3"),
                                          c(ncol(her2n_geo_cg),ncol(her2n_geo2_cg),
                                            ncol(her2n_geo3_cg)))))
rownames(anno_col)=colnames(diff)

pheatmap(diff[sample(1:nrow(diff),1000),],scale="none",annotation_col=anno_col,#annotation_colors=ann_color,
         clustering_method = 'ward.D2',
         show_colnames=F,show_rownames=F,
         fontsize_row=6,fontsize=10,
         legend_breaks=c(0.2,0.8),legend_labels=c('low','high'),
         cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c('#436eee','white','#EE0000'))(100)
)

#去批次
library(dplyr)
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")

#her2p
 
load("F:/002/code/data/her2p/her2p_geo_cg.RData")
load("F:/002/code/data/her2p/her2p_geo2_cg.RData")
load("F:/002/code/data/her2p/her2p_geo3_cg.RData")
load("F:/002/code/data/her2p/her2p_geo4_cg.RData")

cg = intersect(intersect(rownames(her2p_geo_cg),rownames(her2p_geo2_cg)),
                 intersect(rownames(her2p_geo3_cg),rownames(her2p_geo4_cg)))

her2p_geo_cg = her2p_geo_cg[cg,]
her2p_geo2_cg = her2p_geo2_cg[cg,]
her2p_geo3_cg = her2p_geo3_cg[cg,]
her2p_geo4_cg = her2p_geo4_cg[cg,]

methy = cbind(her2p_geo_cg,her2p_geo2_cg,her2p_geo3_cg,her2p_geo4_cg)
sample = colnames(methy)
gse = rep(c("geo1","geo2","geo3","geo4"),
      c(ncol(her2p_geo_cg),ncol(her2p_geo2_cg),
      ncol(her2p_geo3_cg),ncol(her2p_geo4_cg)))
library(sva)
combat_data =ComBat(dat=methy, batch=gse,par.prior=TRUE,prior.plots=FALSE)#校正批次效应

her2p_geoc_cg = as.data.frame(combat_data)
save(her2p_geoc_cg,file = "F:/002/code/data/her2p/her2p_geoc_cg.RData")

her2p_geoc_cg$cg = rownames(her2p_geoc_cg)

her2p_geoc_gene = merge(cgtest,her2p_geoc_cg,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2p_geoc_gene = her2p_geoc_gene[!(grepl(";",her2p_geoc_gene$V3) | her2p_geoc_gene$V3 == "" ),]
her2p_geoc_gene =  as.data.frame(her2p_geoc_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2p_geoc_gene = her2p_geoc_gene[,-(which(colnames(her2p_geoc_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2p_geoc_gene) = her2p_geoc_gene[,1]
her2p_geoc_gene = her2p_geoc_gene[,-1]

save(her2p_geoc_gene,file = "F:/002/code/data/her2p/her2p_geoc_gene.RData")

#her2n

load("F:/002/code/data/her2n/her2n_geo_cg.RData")
load("F:/002/code/data/her2n/her2n_geo2_cg.RData")
load("F:/002/code/data/her2n/her2n_geo3_cg.RData")

cg = intersect(intersect(rownames(her2n_geo_cg),rownames(her2n_geo2_cg)),
                 rownames(her2n_geo3_cg))

her2n_geo_cg = her2n_geo_cg[cg,]
her2n_geo2_cg = her2n_geo2_cg[cg,]
her2n_geo3_cg = her2n_geo3_cg[cg,]

methy = cbind(her2n_geo_cg,her2n_geo2_cg,her2n_geo3_cg)
sample = colnames(methy)
gse = rep(c("geo1","geo2","geo3"),c(ncol(her2n_geo_cg),
                                    ncol(her2n_geo2_cg),ncol(her2n_geo3_cg)))
library(sva)
combat_data =ComBat(dat=methy, batch=gse,par.prior=TRUE,prior.plots=FALSE)#校正批次效应

her2n_geoc_cg = as.data.frame(combat_data)
save(her2n_geoc_cg,file = "F:/002/code/data/her2n/her2n_geoc_cg.RData")

her2n_geoc_cg$cg = rownames(her2n_geoc_cg)
her2n_geoc_gene = merge(cgtest,her2n_geoc_cg,by.x = "V1",by.y = "cg",all.x = F,all.y = F)
her2n_geoc_gene = her2n_geoc_gene[!(grepl(";",her2n_geoc_gene$V3) | her2n_geoc_gene$V3 == "" ),]
her2n_geoc_gene =  as.data.frame(her2n_geoc_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

her2n_geoc_gene = her2n_geoc_gene[,-(which(colnames(her2n_geoc_gene) %in% c("V1","V2","V4","V5")))]
rownames(her2n_geoc_gene) = her2n_geoc_gene[,1]
her2n_geoc_gene = her2n_geoc_gene[,-1]

save(her2n_geoc_gene,file = "F:/002/code/data/her2n/her2n_geoc_gene.RData")
