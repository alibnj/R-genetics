source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
require(GenomicRanges)

#--------------------------------------------------------
#   18s - - - 3UTR
#--------------------------------------------------------

rm183 <- read.delim('C:/Users/banijamali.s/Desktop/N Blast/18S_3UTR_5WS_200mEVAL_iden100-ONLY 7-ONLY 3UTR.txt', header=TRUE, sep = "\t")
rm183 <- rm183[c('Query.Label', 'Target.Label', 'Start.in.Query', 'End.in.Query', 'Start.in.Target', 'E.Value', 'Bit.Score')]
rm183$Query.Label <- '18S'
rm183$Start.in.Target <- '+'
colnames(rm183) <- c('Seq', 'mRNA', 'Start', 'End', 'Strand', 'E.value', 'Bit.Scor')
rm183 <- unique(rm183)


comp183 <- GRanges(seqname=c('18S'), IRanges(start=c(1:1863), end=c(7:1869), names=c(1:1863)))
seqrange183 <- IRanges(start=c(1:1863), end=c(7:1869), names=c(1:1863))


gr183 <- GRanges(seqnames=rm183$Seq, ranges=IRanges(start=rm183$Start, end=rm183$End, names=rm183$mRNA), strand=rm183$Strand)
seq.frq.183 <- as.data.frame(seqrange183)
seq.frq.183$count <- countOverlaps(comp183, gr183, type = "equal")
write.csv(seq.frq.183, 'C:/Users/banijamali.s/Dropbox/Genetics/R/New Results/seq.frq.183-(7)-Eq.csv')


ov183 <- findOverlaps(comp183, gr183, type = "equal")
match_hit <- data.frame(names(comp183)[queryHits(ov183)], names(gr183)[subjectHits(ov183)], stringsAsFactors=F)
names(match_hit) <- c('rRNA','mRNA')
write.csv(match_hit, 'C:/Users/banijamali.s/Dropbox/Genetics/R/New Results/183Match_Hits.csv')


#--------------------------------------------------------
#   18s - - - 5UTR
#--------------------------------------------------------

rm185 <- read.delim('C:/Users/banijamali.s/Desktop/N Blast/18S_5UTR_5WS_120mEVAL_iden100-ONLY 6-ONLY 5UTR.txt', header=TRUE, sep = "\t")
rm185 <- rm185[c("Query.Label", "Target.Label", "Start.in.Query", "End.in.Query", "Start.in.Target", "E.Value", "Bit.Scor") ]
rm185$Query.Label <- '18S'
rm185$Start.in.Target <- '+'
colnames(rm185) <- c('Seq', 'mRNA', 'Start', 'End', 'Strand', 'E.value', 'Bit.Scor')
rm185 <- unique(rm185)


comp185 <- GRanges(seqname=c('18S'), IRanges(start=c(1:1864), end=c(6:1869), names=c(1:1864)))
seqrange185 <- IRanges(start=c(1:1864), end=c(6:1869), names=c(1:1864))


gr185 <- GRanges(seqnames=rm185$Seq, ranges=IRanges(start=rm185$Start, end=rm185$End, names=rm185$mRNA), strand=rm185$Strand)
seq.frq.185 <- as.data.frame(seqrange185)
seq.frq.185$count <- countOverlaps(comp185, gr185, type = "equal")
write.csv(seq.frq.185, 'C:/Users/banijamali.s/Dropbox/Genetics/R/New Results/seq.frq.185-(6)-2Eq.csv')


ov185 <- findOverlaps(comp185, gr185, type = "equal")
match_hit <- data.frame(names(comp185)[queryHits(ov185)], names(gr185)[subjectHits(ov185)], stringsAsFactors=F)
names(match_hit) <- c('rRNA','mRNA')
write.csv(match_hit, 'C:/Users/banijamali.s/Dropbox/Genetics/R/New Results/185Match_Hits.csv')
