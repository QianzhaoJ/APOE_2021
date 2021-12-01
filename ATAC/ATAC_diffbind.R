###############
library(getopt)
command=matrix(c("sample","s",1,"character",'sample prefix',
                 "wt","w",1,"character",'wildtype prefix',
                 "peak","p",1,"character",'peak path',
                 "bam","b",1,"character",'bam path',
                 "number","n",1,"character",'number of replications',
                 "out","o",1,"character",'output path',
                 "help","h",0,"logical",'Help documentation'),
                  byrow=T,ncol=5)
args=getopt(command)
if(!is.null(args$help) ||is.null(args$sample) ||is.null(args$wt) ||is.null(args$peak) || is.null(args$bam) ||is.null(args$out) ||is.null(args$number)) {
  cat(paste(getopt(command, usage = T),"\n"))
  cat('Funtion: diff peak using diffbind',"\n")
  q()
}

###############################################
sample_pre=args$sample
wt_pre=args$wt
peak.path=args$peak
bam.path=args$bam
n=as.numeric(args$number)
out.path=args$out
############################ sample_information
setwd(out.path)
sample_info=data.frame(sample=c(paste(sample_pre,'_rep',1:n,sep = ''),paste(wt_pre,'_rep',1:n,sep = '')),
                       Factor=c(rep(sample_pre,n),rep(wt_pre,n)),
                       Replicate=rep(1:n,2),
                       bamReads=c(paste(bam.path,sample_pre,'_',1:n,'//',sample_pre,'_',1:n,'_no_chrY_chrM.sort.rmdup.bam',sep = ''),paste(bam.path,wt_pre,'_',1:n,'//',wt_pre,'_',1:n,'_no_chrY_chrM.sort.rmdup.bam',sep = '') ),  
                       Peaks=c(paste(peak.path,sample_pre,'_',1:n,'//',sample_pre,'_',1:n,'_peaks.narrowPeak',sep = ''),paste(peak.path,wt_pre,'_',1:n,'//',wt_pre,'_',1:n,'_peaks.narrowPeak',sep = '') ),
                       PeakCaller=rep('narrowPeak',n)
                       
                       )
write.csv(sample_info,'SampleSheet.csv',row.names = F)

#############################
library('data.table')
library('DiffBind')
dbObj <- dba(sampleSheet='SampleSheet.csv')
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

#################
pdf(paste('01_PCA_all_peaks.pdf',sep = ''),width = 8,height = 6)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()
pdf(paste('02_heatmap_all_peaks_correlation.pdf',sep = ''),width = 8,height = 10)
plot(dbObj)
a=plot(dbObj)
write.csv(a,'02_heatmap_all_peaks_correlation.csv',row.names=F)
dev.off()

# Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
# summary of results
re_sum=dba.show(dbObj, bContrasts=T)
write.csv(re_sum,'result_summary.csv')
############################## plot

pdf(paste('03_PCA_sig_peak.pdf',sep = ''),width = 8,height = 6)
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
dev.off()
pdf(paste('04_Venn_comparsion_DEseq2_edgeR.pdf',sep = ''),width = 6,height = 6)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()
pdf(paste('05_PlotMA_',sample_pre,'_vs_',wt_pre,'.pdf',sep = ''),width = 8,height = 6)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
dev.off()
pdf(paste('06_Plobox_',sample_pre,'_vs_',wt_pre,'.pdf',sep = ''),width = 8,height = 6)
dba.plotBox(dbObj)
dev.off()
#####################  write
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file=paste(sample_pre,'_vs_',wt_pre,"_edgeR.txt",sep = '') ,sep="\t", quote=F, col.names = NA)

# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file=paste(sample_pre,'_vs_',wt_pre,"_deseq2.txt",sep = '') , sep="\t", quote=F, col.names = NA)
