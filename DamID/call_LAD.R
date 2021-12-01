############################## APOE LAD
#####################
library(data.table)
library(HMMt)
library(GenomicRanges)
options(scipen = 20)
################ 
out.path="xxxx"
setwd(out.path)
chr.list=paste("chr",c(1:22,"X"),sep = "")
ratio.df=fread("APOE_EMDvsDam_bin25000.bedGraph")
colnames(ratio.df)=c("chr", "start", "end", "score")
ratio.df=ratio.df[ratio.df$chr %in% chr.list,]

ratio.range=GRanges(
  seqnames = ratio.df$chr,
  ranges = IRanges(start = ratio.df$start+1,end = ratio.df$end ),
  score = ratio.df$score
)

######################
bin25.df=fread('G:\\Source\\bin\\bin_25000.bed')
bin25.df=bin25.df[bin25.df$V1 %in% chr.list,]
bin25.range=GRanges(
  seqnames = bin25.df$V1,
  ranges = IRanges(start = bin25.df$V2+1,end = bin25.df$V3 )
)

##########################
overlap.df= findOverlapPairs(bin25.range,ratio.range)@first
overlap.df$score=findOverlapPairs(bin25.range,ratio.range)@second$score
overlap.df=data.frame(overlap.df)[,c(1,2,3,6)]
colnames(overlap.df)=c("chr","start","end","signals")
###########################

bin25.df=bin25.df[,c(1,2,3)]
bin25.df$V2=bin25.df$V2+1
colnames(bin25.df)=c("chr","start","end")

signal.df=merge(bin25.df,overlap.df,by=c("chr","start","end"),all.x=T)
signal.df[is.na(signal.df$signals),]$signals=0
############################

hmm_calls.1 <- HMM(signal.df, na_solution = "NA",file ="APOE.stat" )
hmm_ranges.1 <- getHMMRanges(as(hmm_calls.1, "GRanges"), score = "LAD")
hmm_ranges.2 <- getHMMRanges(as(hmm_calls.1, "GRanges"), score = "iLAD")

write.table(data.frame(hmm_ranges.1)[,1:3],"APOE.HMMT.LAD.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(data.frame(hmm_ranges.2)[,1:3],"APOE.HMMT.iLAD.bed",row.names = F,col.names = F,quote = F,sep = "\t")