######Video source: https://shop119322454.taobao.com
#source("http://bioconductor.org/biocLite.R")
#BiocManager::install(c("DESeq"))
#install.packages("gplots")

foldChange=1
padj=0.05


setwd("/home/kyro_zhang/ZQX/rename")                   #
library("DESeq")
library("limma")

rt=read.table("mRNA.symbol.txt",sep="\t",header=T,check.names=F)  
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>2,] 
data=round(data,0)
nrow(data)

group=c(rep("normal",113),rep("tumor",1109))   #
design = factor(group)
newTab = newCountDataSet( data, design )
newTab = estimateSizeFactors(newTab)
newData=counts(newTab, normalized=TRUE )

#have replicates
newTab = estimateDispersions( newTab, fitType = "local") #估計内部離散度
diff = nbinomTest( newTab, "normal", "tumor") #normal和tumor之間的差異
diff = diff[is.na(diff$padj)==FALSE,]
diff = diff[order(diff$pval),]
write.table( diff, file="DESeqOut_D.csv",sep="\t",quote=F,row.names=F)
diffSig = diff[(diff$padj < padj & (diff$log2FoldChang>foldChange | diff$log2FoldChange<(-foldChange))),]
write.table( diffSig, file="diffSig_D.csv",sep="\t",quote=F,row.names=F)
diffUp = diff[(diff$padj < padj & (diff$log2FoldChang>foldChange)),]
write.table( diffUp, file="up_D.csv",sep="\t",quote=F,row.names=F)
diffDown = diff[(diff$padj < padj & (diff$log2FoldChange<(-foldChange))),]
write.table( diffDown, file="down_D.csv",sep="\t",quote=F,row.names=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp_D.csv",sep="\t",quote=F,col.names=F)   #normalizeExp.txt
diffExp=rbind(id=colnames(newData),newData[diffSig$id,])
write.table(diffExp,file="diffmRNAExp_D.csv",sep="\t",quote=F,col.names=F)         #diffmiRNAExp.txt

#heatmap
hmExp=log10(newData[diffSig$id,]+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",width=60,height=30)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()

#volcano
pdf(file="vol.pdf")
allDiff=diff[is.na(diff$padj)==FALSE,]
max(-log10(allDiff$padj))

xMax=max(-log10(allDiff$padj))+1
yMax=10
plot(-log10(allDiff$padj), allDiff$log2FoldChange, xlab="-log10(padj)",ylab="log2FoldChange",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange>foldChange,]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange<(-foldChange),]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()
