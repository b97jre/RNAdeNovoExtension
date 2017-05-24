library("plyr")
library("ggplot2")
setwd("/Users/johanreimegard/Vetenskap/Data/deNovoExtension/13_2")

sample = "acr_oases" 

initial <- read.table(paste(sample, "initial.info", sep="_"),header=TRUE,sep="\t")
initial$Extension = "Initial"
extended <- read.table(paste(sample, "extended.info", sep="_"),header=TRUE,sep="\t")
extended$Extension = "Extended"
final <- read.table(paste(sample, "final.info", sep="_"),header=TRUE,sep="\t")
final$Extension = "Final"


All =  rbind(initial, extended,final)
All.long = All[All$SeqLength > 350 & All$ORFlength>300, ]
Length  = split(All.long$ORFlength, All.long$Extension)

print("Number of long ORFs")
for(i in 1:length(names(Length))){
  print(paste(names(Length)[i] , length(as.integer(unlist(Length[names(Length)[[i]]]))), sep = " = "))
}

All.ORF = All[All$kind == "FULL" & All$ORFlength > 300 , ]
All.ORF = All[]

ggplot(All.ORF, aes(ORFlength, colour = Extension)) + stat_ecdf()+coord_cartesian(xlim = c(0,5000))

All.all <- summarize(All, ORFlength = unique(ORFlength), 
                           ecdf = ecdf(ORFlength)(unique(ORFlength))* length(ORFlength))

All.Progress <- ddply(All, .(Extension), summarize,
                           ORFlength = unique(ORFlength), 
                           ecdf = ecdf(ORFlength)(unique(ORFlength))* length(ORFlength))

All.ORF.Progress <- ddply(All.ORF, .(Extension), summarize,
                      ORFlength = unique(ORFlength), 
                      ecdf = ecdf(ORFlength)(unique(ORFlength))* length(ORFlength))

All.ORF.Progress$Norm.Ecdf = All.ORF.Progress$ec
head()

ggplot(All.ORF.Progress, aes(ORFlength, ecdf, color = Extension)) + geom_step()

densityPlot = ggplot(All_Long.progress, aes(ORFlength, ..density.., colour = Extension)) 
densityPlot + geom_freqpoly(binwidth = 50)+xlim(500,2000)



diff <- read.table(paste(sample, "initial.round3.diff.ORFs.info", sep="."),header=FALSE,sep="\t")
colnames(diff) <- c("Name","dir","ORFlength","SeqLength","start","stop","kind")
diff$Name <- gsub(">","",diff$Name)
diff$Name <- gsub(" ","",diff$Name)
 diff$kind2 = "ORF"
 diff$kind2[diff$kind==1] = "5UTR"
 diff$kind2[diff$kind==2] = "3UTR"
 diff$kind2[diff$kind==3] = "FULL"
 diff$Extension = "Extended"
 



after <- read.table(paste(sample, "round3.ORFs.info", sep="."),header=FALSE,sep="\t")
colnames(after) <- c("Name","dir","ORFlength","SeqLength","start","stop","kind")
after$Name <- gsub(">","",after$Name)
after$Name <- gsub(" ","",after$Name)
 after$kind2 = "ORF"
after$kind2[after$kind==1] = "5UTR"
after$kind2[after$kind==2] = "3UTR"
after$kind2[after$kind==3] = "FULL"
after$Extension = "After"
 
 
 
 


differentNames  <- setdiff(after$Name,before$Name)

lengths <- lapply(differentNames,nchar)
df <- ldply(lengths)
after$beforeName <- substring(differentNames, first=0, last=df$V1-2)


differentNames1  <- setdiff(beforeName,before$Name)
length(differentNames)

lengths <- lapply(differentNames1,nchar)
df <- ldply(lengths)
beforeName1 <- substring(differentNames1, first=0, last=df$V1-2)

differentNames2  <- setdiff(beforeName1,before$Name)



  
  
merged <- read.table("MergedOnly.info",sep="\t", header=TRUE)
colnames(merged) <- c("Name","dir","ORFlength","SeqLength","start","stop","kind")
merged$Name <- gsub(">","",merged$Name)
merged$Name <- gsub(" ","",merged$Name)
test <- strsplit(x=merged$Name, split="_Merged_")
test2 <- ldply(test)
colnames(test2) <- c("Merged1", "Merged2")
merged$Merged1 <- test2$Merged1
merged$Merged2 <- test2$Merged2

mergedCompare<- merge(merged, after,by.x = "Merged1", by.y ="Name" )
mergedCompare<- merge(mergedCompare, after,by.x = "Merged2", by.y ="Name" )





GCContentBefore <- read.table(paste(sample, "initial.fa.info", sep="."),header=FALSE,sep="\t")
colnames(GCContentBefore) <- c("Name","SeqLength","GC")

GCContentDiff <- read.table(paste(sample, "initial.round3.diff.fa.info", sep="."),header=FALSE,sep="\t")
colnames(GCContentDiff) <- c("Name","SeqLength","GC")

  
  
############################  
#pfam information  START
############################
  
pfamContentBefore <- read.table(paste(sample,"initial_merged.pfam", sep = "."),comment.char="#",  header=FALSE,strip.white=TRUE)
colnames(pfamContentBefore) <- c("seq_id","alignment_start","alignment_end","envelope_start","envelope_end","hmm_acc","hmm_name","type","hmm_start","hmm_end","hmm_length","bit_score","E_value","significance","clan")
pfamContentAfter <- read.table(paste(sample,"round3_merged.pfam", sep = "."),comment.char="#",  header=FALSE,strip.white=TRUE)
colnames(pfamContentAfter) <- c("seq_id","alignment_start","alignment_end","envelope_start","envelope_end","hmm_acc","hmm_name","type","hmm_start","hmm_end","hmm_length","bit_score","E_value","significance","clan")
pfamContentMerged <- read.table("MergedOnly_merged.pfam",comment.char="#",  header=FALSE,strip.white=TRUE)
colnames(pfamContentMerged) <- c("seq_id","alignment_start","alignment_end","envelope_start","envelope_end","hmm_acc","hmm_name","type","hmm_start","hmm_end","hmm_length","bit_score","E_value","significance","clan")

pfamContentBefore$seq_id_Hmm_id <- paste(pfamContentBefore$seq_id,pfamContentBefore$hmm_acc,sep="_")
pfamContentAfter$seq_id_Hmm_id <- paste(pfamContentAfter$seq_id,pfamContentAfter$hmm_acc,sep="_")
pfamContentMerged$seq_id_Hmm_id <- paste(pfamContentMerged$seq_id,pfamContentMerged$hmm_acc,sep="_")
pfamContentMerged$seq_id <- gsub(">","",pfamContentMerged$seq_id)

differentCombos <- setdiff(pfamContentAfter$seq_id_Hmm_id,pfamContentBefore$seq_id_Hmm_id)

pfamContentDiff <- pfamContentAfter[pfamContentAfter$seq_id_Hmm_id %in% differentCombos, ]

#differentNames  <- setdiff(pfamContentAfter$hmm_name,pfamContentBefore$hmm_name)

#pfamContentDiff <- pfamContentAfter[pfamContentAfter$seq_id %in% differentNames, ]

pfamFull <- merge(pfamContentDiff, after,by.x = "seq_id", by.y ="Name")

SameNames <- intersect(pfamContentDiff$seq_id,before$Name)
length(SameNames)



pfamFullCompare<- merge(pfamFull, before,by.x = "seq_id", by.y ="Name" )

pfamFullCompareCOUNT <- merge(pfamFullCompare, idxstats,by.x = "seq_id", by.y ="Name" )

pfamContentAfter()

################################
# PFAM content END
###############################




  
##################################
#BLAST START
#########################################
  initialBlast <- read.table("spruce.13_2.trinity.initial.Picea_Sitchensis.peptide.fa.blast.bestCoverage",header =TRUE,sep ="\t")
  round3Blast <- read.table("spruce.13_2.trinity.round3.Picea_Sitchensis.peptide.fa.blast.bestCoverage",header =TRUE,sep ="\t")
  diffBlast <- read.table("spruce.13_2.trinity.initial.round3.diff.Picea_Sitchensis.peptide.fa.blast.bestCoverage",header =TRUE,sep ="\t")
  
  HitsDatabase <- read.table("Picea_Sitchensis.peptide.fa.size", sep="\t")
  colnames(HitsDatabase) <- c("Name","length")
  
  round3Blast <- round3Blast[round3Blast$BitScore >0 , ]
  
  round3BlastHitsInfo <- merge(round3Blast,HitsDatabase,by.x="Hit", by.y="Name")
  round3BlastHitsInfo <- round3BlastHitsInfo[with(round3BlastHitsInfo,order(Hit,-length.x)), ]
  round3BlastHitsInfo$duplicated <- duplicated(round3BlastHitsInfo$Hit)
  round3BlastHitsInfo <- round3BlastHitsInfo[round3BlastHitsInfo$duplicated==FALSE,  ]
  round3BlastHitsInfo$coverage <- round3BlastHitsInfo$length.x/round3BlastHitsInfo$length.y

  
  
  
  
  initialBlastHits <- initialBlast[initialBlast$BitScore >0 , ]
  
  InitialBlastHitsInfo <- merge(initialBlast,HitsDatabase,by.x="Hit", by.y="Name")
  InitialBlastHitsInfo <- InitialBlastHitsInfo[with(InitialBlastHitsInfo,order(Hit,-length.x)), ]
  InitialBlastHitsInfo$duplicated <- duplicated(InitialBlastHitsInfo$Hit)
  InitialBlastHitsInfo <- InitialBlastHitsInfo[InitialBlastHitsInfo$duplicated==FALSE,  ]
  InitialBlastHitsInfo$coverage <- InitialBlastHitsInfo$length.x/InitialBlastHitsInfo$length.y
  
  
  round3BlastHitsInfoSS <- round3BlastHitsInfo[round3BlastHitsInfo$Hit %in% InitialBlastHitsInfo$Hit , ]
  InitialBlastHitsInfoSS <- InitialBlastHitsInfo[InitialBlastHitsInfo$Hit %in% round3BlastHitsInfoSS$Hit , ]
  
  plot(InitialBlastHitsInfoSS$length.x,round3BlastHitsInfoSS$length.x)
  histInitial <- hist(InitialBlastHitsInfo$coverage)
  histRond3 <- hist(round3BlastHitsInfo$coverage)
   
 
  histInitial <- hist(InitialBlastHitsInfo$coverage)
  histRond3 <- hist(round3BlastHitsInfo$coverage)
  
 heading = "Best hit coverage on cDNA library"
 plot(histInitial$mids,histInitial$counts, type="n", main=heading)
 lines(histInitial$mids,histInitial$counts)
 lines(histRond3$mids,histRond3$counts, col = "blue")
 
 
 
  plot(density(diffBlastInfo$coverage))
lines(density((InitialBlastHitsInfosubset$coverage)))
 plot(InitialBlastHitsInfosubset$length.x,diffBlastInfo$length.x)
 
  
  
  plot(diffBlastInfoBestHit$length.x,diffBlastInfoBestHit$length.y)
  hist(diffBlastInfoBestHit$coverage)
  diffBlastHits <- diffBlast[diffBlast$BitScore > 0, ] 
  initialBlastHitsDiff <- initialBlastHits[initialBlastHits$Query  %in% diffBlastHits$Query, ]
  diffBlastHitsInitial  <- diffBlastHits[diffBlastHits$Query %in% initialBlastHitsDiff$Query, ]
 
  beforeGeneInfo <- before[before$Name %in%  initialBlastHitsDiff$Query, ]
  afterGeneInfo <- after[after$Name %in%  initialBlastHitsDiff$Query, ]
  
  # sort transcripts
 
  diffBlastHitsInitial <- diffBlastHitsInitial[with(diffBlastHitsInitial,order(Query)), ]
  initialBlastHitsDiff <- initialBlastHitsDiff[with(initialBlastHitsDiff,order(Query)), ]
 
 beforeGeneInfo <- beforeGeneInfo[with(beforeGeneInfo,order(Name)), ]
  afterGeneInfo <- afterGeneInfo[with(afterGeneInfo,order(Name)), ]
  
  
  
  afterGeneInfo$diff <- afterGeneInfo$ORFlength - beforeGeneInfo$ORFlength
  GeneDiff <-   afterGeneInfo[afterGeneInfo$diff > 100 , ]
 
  beforeGeneInfoDiff <- beforeGeneInfo[beforeGeneInfo$Name %in% GeneDiff$Name, ]
  initialBlastHitsDiffOpen <- initialBlastHitsDiff[initialBlastHitsDiff$Query  %in% GeneDiff$Name, ]
  diffBlastHitsInitialOpen  <- diffBlastHitsInitial[diffBlastHitsInitial$Query %in% GeneDiff$Name, ]
  
  
  
  diffBlastHitsInitialOpen$HitDiff<- c(diffBlastHitsInitialOpen$length - initialBlastHitsDiffOpen$length)
  diffBlastHitsInitialOpen$ORFDiff <- c((GeneDiff$ORFlength - beforeGeneInfoDiff$ORFlength)/3)
  diffBlastHitsInitialOpen$Percentage <-   diffBlastHitsInitialOpen$HitDiff/diffBlastHitsInitialOpen$ORFDiff
  
 plot(density)
  
  diffBlastHitsInitialOpenStrange <- diffBlastHitsInitialOpen[diffBlastHitsInitialOpen$Percentage<1.5, ]
  HitsDatabaseStrange <- HitsDatabase[HitsDatabase$Name %in% diffBlastHitsInitialOpenStrange$Hit, ]
  
  
  
  
  
  
  
  
  head()
  
  par(mfrow=c(2,1), bty="l", cex=0.6)
  plot(initialBlastHitsDiffOpen$BitScore,diffBlastHitsInitialOpen$BitScore)
  plot(diffBlastHitsInitialOpen$HitDiff,diffBlastHitsInitialOpen$ORFDiff)
  
  
  
  
########################################
#BLAST STOP
#######################################
  
#########################################
# Merging info from different sources START
##########################################


transcripts <- union(mergedCompare$Merged2, mergedCompare$Merged1)

pfamContentDiffMerged <- pfamContentDiff[pfamContentDiff$seq_id %in% transcripts, ]

pfamContentDiffMerged1 <- merge(pfamContentDiffMerged, mergedCompare, by.x="seq_id" , by.y= "Merged1")
pfamContentDiffMerged2 <- merge(pfamContentDiffMerged, mergedCompare, by.x="seq_id" , by.y= "Merged2")

pfamContentDiffMerged1Total <- merge(pfamContentDiffMerged1, pfamContentMerged, by.x=c("Name","hmm_acc","hmm_name") , by.y= c("seq_id","hmm_acc","hmm_name" ))
pfamContentDiffMerged2Total <- merge(pfamContentDiffMerged2, pfamContentMerged, by.x=c("Name","hmm_acc","hmm_name") ,by.y= c("seq_id","hmm_acc","hmm_name" ))

pfamContentDiffMerged1Total$BitScoreChange = pfamContentDiffMerged1Total$bit_score.y-pfamContentDiffMerged1Total$bit_score.x
pfamContentDiffMerged1Total <- pfamContentDiffMerged1Total[with(pfamContentDiffMerged1Total,order(-BitScoreChange)), ]

pfamContentDiffMerged2Total$BitScoreChange = pfamContentDiffMerged2Total$bit_score.y-pfamContentDiffMerged2Total$bit_score.x
pfamContentDiffMerged2Total <- pfamContentDiffMerged2Total[with(pfamContentDiffMerged2Total,order(-BitScoreChange)), ]


## Should contain pfamDomain 
## Filter 2 SeqLength > 500 & ORFlength > 300
Filtered1 = pfamContentDiffMerged1Total[pfamContentDiffMerged1Total$SeqLength.x >500 & pfamContentDiffMerged1Total$ORFlength.x >300, ]
Filtered2 = pfamContentDiffMerged2Total[pfamContentDiffMerged2Total$SeqLength.x >500 & pfamContentDiffMerged2Total$ORFlength.x >300, ]

## Filter 2 ORFlength > 300
Filtered1_1 = Filtered1[Filtered1$ORFlength.x-Filtered1$ORFlength.y > 100 & Filtered1$ORFlength.x-Filtered1$ORFlength > 100, ]
Filtered2_1 =Filtered1[Filtered2$ORFlength.x-Filtered2$ORFlength.y > 100 & Filtered2$ORFlength.x-Filtered2$ORFlength > 100, ]


##
Filtered1_1_1 = Filtered1_1[(Filtered1_1$hmm_end.y -Filtered1_1$hmm_start.y )/Filtered1_1$hmm_length.y  > 0.8, ]
Filtered1[Filtered1$ORFlength.x-Filtered1$ORFlength.y > 100 & Filtered1$ORFlength.x-Filtered1$ORFlength > 100, ]

FinalTable = data.frame(Filtered1_1_1$Name,Filtered1_1_1$SeqLength.x,Filtered1_1_1$ORFlength.x,Filtered1_1_1$ORFlength, Filtered1_1_1$ORFlength.y,Filtered1_1_1$hmm_name,Filtered1_1_1$hmm_acc,Filtered1_1_1$alignment_start.y,Filtered1_1_1$alignment_end.y,Filtered1_1_1$E_value.y)
colnames(FinalTable) <- c("Name","Seq_length","ORF_length","ORF_length1","ORF_length2","PFAM","PFAM_ID","PFAM_start","PFAM_stop","PFAM_Evalue")

FinalTable2 = data.frame(Filtered2_1$Name,Filtered2_1$SeqLength.x,Filtered2_1$ORFlength.x,Filtered2_1$ORFlength, Filtered2_1$ORFlength.y,Filtered2_1$hmm_name,Filtered2_1$hmm_acc,Filtered2_1$alignment_start.y,Filtered2_1$alignment_end.y,Filtered2_1$E_value.y)
colnames(FinalTable2) <- c("Name","Seq_length","ORF_length","ORF_length1","ORF_length2","PFAM","PFAM_ID","PFAM_start","PFAM_stop","PFAM_Evalue")



FFinalTable[FinalTable$PFAM=="Frigida",]

FinalTable = data.frame(Filtered1_1_1$Name)


write.table()

Filtered1_1_1_
head(pfamContentDiffMerged1Total)


pfamContentDiff <- pfamContentDiff[pfamContentDiff$seq_id %in% transcripts, ]



head(pfamContentDiffMerged1)










#########################################
# Merging info from different sources STOP
##########################################

#########################################
# Different graphs START
##########################################



meanGCcontent <- mean(GCContentDiff$GC )
SDGCcontent <-  sd(GCContentDiff$GC)
GCContentBeforeDiff <- GCContentBefore[GCContentBefore$Name %in% GCContentDiff$Name, ]
GCContentBeforeORFs <- GCContentBefore[GCContentBefore$Name %in% beforeORFsLong$Name, ]

head(idxstats)



pdf("GC_content.pdf")
plot(density(GCContentBefore$GC),xlab="Percent GC",main="GC content distribution")
lines(density(GCContentBeforeDiff$GC),col="blue")
lines(density(GCContentDiff$GC),col="red")
lines(density(GCContentBeforeORFs$GC),col="cyan")

legend("topright",legend=c("All","before","after","long orf"),col=c("black","blue","red","cyan"),lty=1)
dev.off()






CountBefore <- read.table("13_21fastq_13_2.trinity.round0.strict.sam.sorted.idxstats",header=FALSE,sep="\t")
colnames(CountBefore) <- c("Name","SeqLength","Mapped", "NotMapped")
head(CountBefore)
CountBefore$RBBP=CountBefore$Mapped/CountBefore$SeqLength*100

CountAfter <- read.table("spruce_13_2_trinity_initial_spruce.13_2.trinity.round3.strict.sorted.idxstats",header=FALSE,sep="\t")
colnames(CountAfter) <- c("Name","SeqLength","Mapped", "NotMapped")
CountAfter$RBBP=CountAfter$Mapped/CountAfter$SeqLength*100



beforeORFsShort <- before[which(before$kind==0 & before$SeqLength>300 & before$SeqLength<1000), ]
beforeORFsLong <- before[which(before$kind==0 & before$SeqLength>999), ]
beforeNonORFsShort <- before[which(before$kind==3 & before$SeqLength>300 & before$SeqLength<1000), ]
beforeFullLong <- before[which(before$kind==3 & before$ORFlength>1000), ]

beforeORFsLong <- before[which(before$kind==0 & before$SeqLength>999), ]

CountBeforeDiff <- CountBefore[CountBefore$Name %in% GCContentDiff$Name, ]
CountAfterDiff <- CountAfter[CountAfter$Name %in% GCContentDiff$Name, ]

beforeORFsLong <- before[which(before$kind==0 & before$SeqLength>1000), ]
CountAfterORFLong <- CountAfter[CountAfter$Name %in% beforeORFsLong$Name, ]
CountAfterORFshort <- CountAfter[CountAfter$Name %in% beforeORFsShort$Name, ]
CountAfterFullShort <- CountAfter[CountAfter$Name %in% beforeNonORFsShort$Name, ]
CountAfterFullLong <- CountAfter[CountAfter$Name %in% beforeFullLong$Name, ]
CountAfterDiffLong <- CountAfter[CountAfter$Name %in% GCContentDiff$Name, ]

CountAfterDiff <- CountAfter[CountAfter$Name %in% GCContentDiff$Name, ]

ReadCoverageDiffSelected <- CountAfterDiff[ CountAfterDiff$RBBP>25 ,]
ORFlengthDiffSelected <- diff[ diff$ORFlength >1000, ]
GCContentDiffSelected <- GCContentDiff[ GCContentDiff$GC > meanGCcontent-SDGCcontent & GCContentDiff$GC < meanGCcontent+SDGCcontent , ] 

FullfillAll = intersect(ReadCoverageDiffSelected$Name,intersect(ORFlengthDiffSelected$Name,GCContentDiffSelected$Name))

write.table(FullfillAll,file="AllFullilled.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

diffAllFullfilled  <- diff[ diff$Name  %in% FullfillAll,  ]
GCAllFullfilled  <- GCContentDiff[ GCContentDiff$Name  %in% FullfillAll,  ]
ReadCoverageDiffSelected <- CountAfterDiff[ CountAfterDiff$Name  %in% FullfillAll ,]


ReadCoverage = list(ReadCoverage = ReadCoverageDiffSelected$Name, ORFlength = ORFlengthDiffSelected$Name, GCContent = GCContentDiffSelected$Name )
grid(venn.diagram(ReadCoverage,NULL))
ORFlength = list()
GCContent = list(GCContent = GCContentDiffSelected)






pdf("Count_content.pdf")
plot(density(CountAfter$RBBP[complete.cases(CountAfter$RBBP) & CountAfter$RBBP < 100]),xlab="#Reads*100/length",main="Reads density", xlim=c(0,60))
lines(density(CountAfterORFLong$RBBP), col="cyan")
lines(density(CountBeforeDiff$RBBP[CountBeforeDiff$RBBP < 1000] ),col="blue",)
lines(density(CountAfterDiff$RBBP[CountAfterDiff$RBBP < 1000 ]), col="red")
legend("topright",legend=c("All","before","after","long orf"),col=c("black","blue","red","cyan"),lty=1)
dev.off()
 
 lines(density(CountAfterORFLong$RBBP[CountAfterORFLong$RBBP < 100]), col="cyan")
 lines(density(CountAfterORFshort$RBBP[CountAfterORFshort$RBBP < 100]), col="green")
lines(density(CountAfterFullShort$RBBP[CountAfterFullShort$RBBP < 100]), col="yellow")
lines(density(CountAfterFullLong$RBBP[CountAfterFullLong$RBBP < 100]), col="brown")
 
legend("topright",legend=c("All","Long orfs","Extended"),col=c("black","blue","red","cyan"),lty=1)

dev.off()



plot(density(CountAfter$RBBP[complete.cases(CountAfter$RBBP) & CountAfter$RBBP < 50]))
lines(density(CountAfterDiff$RBBP[CountAfterDiff$RBBP < 50]))



pdf("SequenceSize.pdf")

plot(density(idxstats$length),col="black",xlim=c(0,2000),xlab="Sequence length",main="Sequence lengths")
lines(density(before$SeqLength),col="blue")
lines(density(after$SeqLength),col="red")
legend("topright",legend=c("All","before","after"),col=c("black","blue","red"),lty=1)

dev.off()

pdf("ORFSize.pdf")

plot(density(before$ORFlength),col="blue",xlim=c(0,1000), xlab="ORF length",main="ORF lengths")
lines(density(after$ORFlength),col="red")
legend("topright",legend=c("before","after"),col=c("blue","red"),lty=1)

dev.off()

pdf("size_Difference.pdf")
beforeORFs <- before[which(before$kind==0), ]
afterSame <- after[match(before$Name, after$Name), ]
afterORFs <- after[match(beforeORFs$Name, after$Name), ]
changed <- afterORFs[which(afterORFs$kind==3), ]


plot(density(afterSame$ORFlength-before$ORFlength),col="blue",xlim=c(0,500),xlab="Difference in length",main="Difference in length")
lines(density(afterSame$SeqLength-before$SeqLength),col="black")
lines(density(afterORFs$ORFlength-beforeORFs$ORFlength),col="red")
legend("topright",legend=c("Sequence","ORF all","ORFs CDS"),col=c("black","blue","red"),lty=1)




dev.off()

pdf("Kind.pdf")
par(mfrow=c(2,1), bty="l", cex=0.6)
hist(before$kind,freq = FALSE,main="Before")
hist(afterSame$kind,freq = FALSE,main="After")
dev.off()
#########################################
# Different graphs STOP
##########################################

