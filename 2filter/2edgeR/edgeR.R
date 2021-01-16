#biocLite("edgeR")
#install.packages("gplots")

foldChange=1  #fold Change基因差異倍數，增大嚴格
padj=0.05 #P-value adjust減少嚴格

# setwd("D:\\2020\\資料科學\\final-project\\tcga_files\\rename\\edgeR")                   
setwd("/home/kyro_zhang/ZQX/rename")
library("edgeR")
rt=read.table("mRNA.symbol.txt",sep="\t",header=T,check.names=F)  #
rt=as.matrix(rt)
rownames(rt)=rt[,1]   #第一列基因名為rowname
exp=rt[,2:ncol(rt)]  #第二列到最後一列為表達量
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

#group=c("normal","tumor","tumor","normal","tumor")
group=c(rep("normal",113),rep("tumor",1109))                    
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)  #轉化為edgeR對象格式
y <- calcNormFactors(y)       #校正因子
y <- estimateCommonDisp(y)    #估計normal的變異（先估計内部差異程度，再看他們之間的差異是否大於内部差異，如果更大，達到一定水平就可以篩選出來作爲我們所需結果）
y <- estimateTagwiseDisp(y)   #估計tumor的變異
et <- exactTest(y,pair = c("normal","tumor"))     #檢驗
topTags(et)
ordered_tags <- topTags(et, n=100000)  #顯示前10w，篩選之後結果小於10w，即所有基因都顯示

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)  #去除低表達之後的基因
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),] 
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)   #顯著差異的基因（fdr，foldchange篩選后的）
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]     #foldchange大於1上調，小於1下調
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #normalized后的高表達基因
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #normalized后的差異基因

heatmapData <- newData[rownames(diffSig),]

#volcano
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

#heatmap
hmExp=log10(heatmapData+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",width=60,height=90)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()
