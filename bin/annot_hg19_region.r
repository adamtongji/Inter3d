library("ChIPseeker")
library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

args=commandArgs(T)

ctcf1_prefix<-args[1]
ctcf2_prefix<-args[2]
share_peak_prefix<-args[3]

peak = readPeakFile(paste0(ctcf1_prefix,".bed"))
peak_annots = annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=txdb,annoDb = "org.Hs.eg.db")
peak_anno_tab<- as.data.frame(peak_annots@anno)
write.table(peak_anno_tab,file = paste0(ctcf1_prefix,"_annot.txt"),sep = "\t",row.names = F,col.names = T,quote=F)

peak = readPeakFile(paste0(ctcf2_prefix,".bed"))
peak_annots = annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=txdb,annoDb = "org.Hs.eg.db")
peak_anno_tab<- as.data.frame(peak_annots@anno)
write.table(peak_anno_tab,file = paste0(ctcf2_prefix,"_annot.txt"),sep = "\t",row.names = F,col.names = T,quote=F)

peak = readPeakFile(paste0(share_peak_prefix,".bed"))
peak_annots = annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=txdb,annoDb = "org.Hs.eg.db")
peak_anno_tab<- as.data.frame(peak_annots@anno)
write.table(peak_anno_tab,file = paste0(share_peak_prefix,"_annot.txt"),sep = "\t",row.names = F,col.names = T,quote=F)

