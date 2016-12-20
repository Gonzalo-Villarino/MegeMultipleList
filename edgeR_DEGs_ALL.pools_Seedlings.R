rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_JasonReed/HtSeq//")
getwd()
library("edgeR")




####################################################################################################################################
# EDGER (using htseq-counts previously)
# edgeR is based on Gener.Linear Model (GML) - assuming that read counts are distributed according to the negative binomial distribution
####################################################################################################################################

# Import (to select prot. coding genes only) 
gene_anno = read.table("TAIR10.gene.list.txt",sep="\t")
names(gene_anno) = c("id","anno")

# Import and merge the 14 htseq-count tables to start Diff.Expr.Analysis
htseqAllCount_BL = data.frame()

filelist = read.table("file.list.txt",header=F)
for(i in 1:nrow(filelist)){
  data = read.table(as.character(filelist[i,1]),sep="\t",header=F)
  names(data)[1] = "id"
  names(data)[2] = strsplit(as.character(filelist[i,1]),"_")[[1]][1]
  if(i==1){
    htseqAllCount_B = data
  } else{
    htseqAllCount_B = merge(htseqAllCount_B,data,by=c("id"))
  }
  
}




##########################################################################

###get protein coding genes###
htseqAllCount_BL <- merge(htseqAllCount_B,gene_anno,by=c(1))
htseqAllCount_BL <- htseqAllCount_BL[htseqAllCount_BL$anno=="protein_coding_gene",]
htseqAllCount_BL = htseqAllCount_BL[,-29]
####


# make table of only 6 numeric values for downstream analysis
cm <- htseqAllCount_BL[,-1]
rownames(cm) <- htseqAllCount_BL[,1]

# build DGEList
group <- c(1,1,1,1,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,7,7,8,8,8,8)
y <- DGEList(counts = cm, group=c(1,1,1,1,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,7,7,8,8,8,8))
str(y)
dim(y)


# paramenter for filtering low expressed genes
min.cpm <- 2
n.min.cpm <- 3
keep <- rowSums(cpm(y)>min.cpm) >= n.min.cpm
table(keep)
y <- y[keep,]
dim(y) #19,233 genes expressed =  27 samples 

y$samples$lib.size <- colSums(y$counts)


# Check distribution
#min.cpm
#head(cpm(y))
#hist(cpm(y)[,1])
#hist(cpm(y)[,1],breaks = 20)
#hist(log(cpm(y)[,1]),breaks = 20)
#hist(log(cpm(y)[,1])/log10,breaks = 20)
#hist(log(cpm(y)[,1])/log(10),breaks = 20)


# TMM normalization
y <- calcNormFactors(y, method="TMM")
y$samples

# prepare for edgeR glm
PRT  <- factor(c("poolA","poolA","poolA","poolA","poolB",
                 "poolB","poolC","poolC","poolC","poolC", "poolD",
                 "poolD","poolD","poolD","poolE","poolE","poolE","poolE","poolF","poolF", 
                 "poolF", "poolG","poolG","poolH","poolH","poolH",
                 "poolH"), 
               levels= c("poolA", "poolB","poolC", "poolD","poolE", "poolF","poolG", "poolH"))
sample.names <- c("poolA_0","poolA_1","poolA_2","poolA_3","poolB_0",
                  "poolB_1","poolC_0","poolC_1","poolC_2","poolC_3", "poolD_0",
                  "poolD_1","poolD_2","poolD_3","poolE_0","poolE_1", "poolE_2","poolE_3","poolF_0","poolF_1", 
                  "poolF_2", "poolG_0","poolG_1","poolH_0","poolH_1","poolH_2","poolH_3")

targets <- as.data.frame(cbind(sample.names,PRT))
design <- model.matrix(~0+PRT)


my.contrasts <- makeContrasts(
        AB = (PRTpoolA-PRTpoolB), # comparison: yfp pos / yfp neg
        AC = (PRTpoolA-PRTpoolC),  # comparison all sort / no sort 
        AD = (PRTpoolA-PRTpoolD), # comparison all sort / yfp pos
        AE = (PRTpoolA-PRTpoolE), # comparison all sort / yfp neg
        DE = (PRTpoolD-PRTpoolE),
        FG = (PRTpoolF-PRTpoolG),
        FH = (PRTpoolF-PRTpoolH), levels=design # comparison no sort / yfp neg
)
interesting.contrasts <- c("AB", "AC", "AD", "AE", "DE", "FG","FH")


# variance estimate
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)


# QC plots
plotMDS(y) # MDS plot
plotMeanVar(y) # Mean vs Variance
plotBCV(y) #


# EdgeR analysis & report
fit <- glmFit(y,design)
fit

# fdr threshold
fdr.t <- 1
project.name <- "RNAseq_edgeR_JReed"
cat("experiment: ", project.name, "\n")
cat("thresholds: ", min.cpm, " ", n.min.cpm,"\n")
for (my.contrast in interesting.contrasts) {
        lrt <- glmLRT(fit, contrast=my.contrasts[,my.contrast])
        etable <- topTags(lrt, n=nrow(lrt$table), adjust.method="BH")
        etable <- etable$table[etable$table$FDR<fdr.t,]
        etable <- etable[ order(etable$FDR), ]  
        cat(my.contrast," ", dim(etable)[1], " genes\n")
        # write result (c1,...,c6) in workign dir. 
        write.table( etable[,], file=paste(project.name,my.contrast, min.cpm,n.min.cpm,sep="."), row.names=TRUE)
}

#### ENDS HERE ####

et <- exactTest(y)
topTags(et)


detags <- rownames(topTags(et)$table)
cpm(y)[detags, order(y$samples$group)]
summary(de <- decideTestsDGE(et, p=0.01))

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "blue")


################################################
#RPKM CALCULTING ( TAIR10.gene.length) #### 
################################################

# calculate gene lenght 
len = read.table("tair10.whole.genelength.txt",sep="\t")
names(len) = c("ID","length")
htseqAllCount_BL = merge(htseqAllCount_BL,len,by=c(1))
x = htseqAllCount_BL[,c(2:29)] # 
rownames(x) = htseqAllCount_BL[,1]

rpkm = rpkm(x,htseqAllCount_BL[,29],normalized.lib.sizes=TRUE)
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

str(rpkm)

### Merge Table ### 

# merge p-values/FDR table with RPKM table
etable$name = rownames(etable)
edgeR_table = merge(rpkm,etable,by=c("name"))
str(table)


#final table
edgeR_table 
write.table(edgeR_table, file="edgeRTable_allSamples_JR")

