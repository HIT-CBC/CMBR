#Boruta
#her2p
library(Boruta)

load("F:/002/code/data/her2p/methy_her2p.RData")
load("F:/002/code/data/methy_normal.RData")
brown = read.table(file = "F:/002/code/data/brown_gene.txt",sep = "\n")
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")

train = sample(methy_her2p,round(ncol(methy_her2p)*0.7),replace = F)
test = methy_her2p[,which(colnames(methy_her2p) %in% setdiff(colnames(methy_her2p),colnames(train)))]

diff = cbind(train,methy_normal)

Pvalue = apply(diff,1,function(x){
t.test(x[1:ncol(train)],x[(ncol(train)+1):ncol(diff)])$p.value
})
gap = apply(train,1,mean)-apply(methy_normal,1,mean)
fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))
train_cg = train[which(fdr <= 0.05 & abs(gap) >=0.2),]

type=rep(c("her2p","normal"),c(ncol(train),ncol(methy_normal)))

methy_normal = cbind(rownames(methy_normal),methy_normal)
colnames(methy_normal)[1] = "cg"

train_cg = cbind(rownames(train_cg),train_cg)
colnames(train_cg)[1] = "cg"

tdifftrain1 = merge(train_cg,methy_normal,by = "cg",all.x = F,all.y = F)
methy_normal = methy_normal[,-1]

rownames(tdifftrain1) = tdifftrain1[,"cg"]
tdifftrain1 = as.data.frame(t(tdifftrain1[,-1]))

tdifftrain1=cbind(tdifftrain1,type)
tdifftrain1[,ncol(tdifftrain1)]=as.factor(tdifftrain1[,ncol(tdifftrain1)])

Boruta(type~.,data=tdifftrain1,doTrace=2)->Boruta.train.extended

boruta_signif <- names(Boruta.train.extended$finalDecision[Boruta.train.extended$finalDecision %in% c("Confirmed", "Tentative")])
boruta_gene = unique(cgtest$V3[which(cgtest$V1 %in% boruta_signif)])

mygene = intersect(boruta_gene,brown[,1])
write.table(mygene,file = "F:/002/code/data/her2p/mygene.txt",quote = F,row.names = F,col.names = F)
save(train,test,boruta_signif,file = "F:/002/code/data/her2p/train_test.RData")

#her2n
library(Boruta)
load("F:/002/code/data/her2n/methy_her2n.RData")
load("F:/002/code/data/methy_normal.RData")
brown = read.table(file = "F:/002/code/data/brown_gene.txt",sep = "\n")
cgtest = read.table(file ="F:/002/code/data/cgtest.txt",header = F,sep = "\t")

train = sample(methy_her2n,round(ncol(methy_her2n)*0.7),replace = F)
test = methy_her2n[,which(colnames(methy_her2n) %in% setdiff(colnames(methy_her2n),colnames(train)))]

diff = cbind(train,methy_normal)

Pvalue = apply(diff,1,function(x){
t.test(x[1:ncol(train)],x[(ncol(train)+1):ncol(diff)])$p.value
})
gap = apply(train,1,mean)-apply(methy_normal,1,mean)
fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))
train_cg = train[which(fdr <= 0.05 & abs(gap) >=0.2),]

type=rep(c("her2n","normal"),c(ncol(train),ncol(methy_normal)))

methy_normal = cbind(rownames(methy_normal),methy_normal)
colnames(methy_normal)[1] = "cg"

train_cg = cbind(rownames(train_cg),train_cg)
colnames(train_cg)[1] = "cg"

tdifftrain1 = merge(train_cg,methy_normal,by = "cg",all.x = F,all.y = F)
methy_normal = methy_normal[,-1]

rownames(tdifftrain1) = tdifftrain1[,"cg"]
tdifftrain1 = as.data.frame(t(tdifftrain1[,-1]))

tdifftrain1=cbind(tdifftrain1,type)
tdifftrain1[,ncol(tdifftrain1)]=as.factor(tdifftrain1[,ncol(tdifftrain1)])

Boruta(type~.,data=tdifftrain1,doTrace=2)->Boruta.train.extended

boruta_signif <- names(Boruta.train.extended$finalDecision[Boruta.train.extended$finalDecision %in% c("Confirmed", "Tentative")])
boruta_gene = unique(cgtest$V3[which(cgtest$V1 %in% boruta_signif)])

mygene = intersect(boruta_gene,brown[,1])

write.table(mygene,file = "F:/002/code/data/her2n/mygene.txt",quote = F,row.names = F,col.names = F)
save(train,test,boruta_signif,file = "F:/002/code/data/her2n/train_test.RData")
######################################################################
#分类器
library(caret)

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

load("F:/002/code/data/methy_cancer.RData")

#her2p
load("F:/002/code/data/her2p/train_test.RData")
load("F:/002/code/data/her2p/methy_her2p.RData")
unher2p = setdiff(colnames(methy_cancer),colnames(methy_her2p))
unher2p_train = sample(unher2p,round(length(unher2p)*0.1),replace = F)
unher2p_test = sample(setdiff(unher2p,unher2p_train),round(length(unher2p)*0.05),replace = F)

traindata = methy_cancer[which(rownames(methy_cancer) %in% boruta_signif),which(colnames(methy_cancer) %in% union(colnames(train),unher2p_train))]
testdata = methy_cancer[which(rownames(methy_cancer) %in% boruta_signif),which(colnames(methy_cancer) %in% union(colnames(test),unher2p_test))]

traindata = as.data.frame(t(traindata))
testdata = as.data.frame(t(testdata))

traindata$type = ifelse(rownames(traindata) %in% unher2p ,"unher2p","her2p")
testdata$type = ifelse(rownames(testdata) %in% unher2p ,"unher2p","her2p")

traindata$type=as.factor(traindata$type)
testdata$type=as.factor(testdata$type)

inTrain<- createDataPartition(y=traindata$type,p=0.75,list=FALSE)
training<- traindata[inTrain, ]
testing<- traindata[-inTrain, ]

modelFit=train(type~.,data=training,method="ranger",metric = "ROC",trControl = fitControl) 

atesting = confusionMatrix(predict(modelFit,newdata=testing),testing$type)$byClass['Balanced Accuracy']
atraining = confusionMatrix(predict(modelFit,newdata=training),training$type)$byClass['Balanced Accuracy']
atest = confusionMatrix(predict(modelFit,newdata=testdata),testdata$type)$byClass['Balanced Accuracy']

#her2n
load("F:/002/code/data/her2n/train_test.RData")
load("F:/002/code/data/her2n/methy_her2n.RData")
unher2n = setdiff(colnames(methy_cancer),colnames(methy_her2n))
unher2n_train = sample(unher2n,round(length(unher2n)*0.7),replace = F)
unher2n_test = setdiff(unher2n,unher2n_train)

traindata = methy_cancer[which(rownames(methy_cancer) %in% boruta_signif),which(colnames(methy_cancer) %in% union(colnames(train),unher2n_train))]
testdata = methy_cancer[which(rownames(methy_cancer) %in% boruta_signif),which(colnames(methy_cancer) %in% union(colnames(test),unher2n_test))]

traindata = as.data.frame(t(traindata))
testdata = as.data.frame(t(testdata))

traindata$type = ifelse(rownames(traindata) %in% unher2n ,"unher2n","her2n")
testdata$type = ifelse(rownames(testdata) %in% unher2n ,"unher2n","her2n")

traindata$type=as.factor(traindata$type)
testdata$type=as.factor(testdata$type)

inTrain<- createDataPartition(y=traindata$type,p=0.75,list=FALSE)
training<- traindata[inTrain, ]
testing<- traindata[-inTrain, ]

modelFit=train(type~.,data=training,method="ranger",metric = "ROC",trControl = fitControl) 

atesting = confusionMatrix(predict(modelFit,newdata=testing),testing$type)$byClass['Balanced Accuracy']
atraining = confusionMatrix(predict(modelFit,newdata=training),training$type)$byClass['Balanced Accuracy']
atest = confusionMatrix(predict(modelFit,newdata=testdata),testdata$type)$byClass['Balanced Accuracy']
