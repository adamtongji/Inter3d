library(stringr)
args=commandArgs(T)

tad1<-args[1]
tad2<-args[2]
outputdir<-args[3]

tad_1_file<-read.table(tad1,sep='\t')
tad_2_file<-read.table(tad2,sep='\t')

tad_1_pos<- unique(c(paste(tad_1_file[,1],tad_1_file[,2],sep="_"),paste(tad_1_file[,1],tad_1_file[,3],sep="_")))
tad_2_Pos<- unique(paste(tad_2_file[,1],tad_2_file[,2],sep="_"),paste(tad_2_file[,1],tad_2_file[,3],sep="_"))

tad_1_spec<-setdiff(tad_1_pos,tad_2_Pos)
tad_2_spec<- setdiff(tad_2_Pos,tad_1_pos)
tad_share <- intersect(tad_2_Pos,tad_1_pos)
# 对于spec分布左右延伸5kb,10kb，看是否有交集(boundary distance)
tad_1_spec_mat <- data.frame(str_split(tad_1_spec,pattern = "_",simplify = T))
tad_1_spec_mat[,2]= as.numeric(as.matrix(tad_1_spec_mat[,2]))
tad_2_spec_mat <- data.frame(str_split(tad_2_spec,pattern = "_",simplify = T))
tad_2_spec_mat[,2]= as.numeric(as.matrix(tad_2_spec_mat[,2]))

# tad_1_spec_5k <-data.frame(chr=tad_1_spec_mat[,1],start=tad_1_spec_mat[,2]-5000,end=tad_1_spec_mat[,2]+5000)
tad_1_spec_10k <-data.frame(chr=tad_1_spec_mat[,1],start=tad_1_spec_mat[,2]-10000,end=tad_1_spec_mat[,2]+10000)
tad_1_spec_10k$start[tad_1_spec_10k$start<1]<- 1
# tad_2_spec_5k <-data.frame(chr=tad_2_spec_mat[,1],start=tad_2_spec_mat[,2]-5000,end=tad_2_spec_mat[,2]+5000)
tad_2_spec_10k <-data.frame(chr=tad_2_spec_mat[,1],start=tad_2_spec_mat[,2]-10000,end=tad_2_spec_mat[,2]+10000)
tad_2_spec_10k$start[tad_2_spec_10k$start<1]<- 1

write.table(tad_1_spec_10k,file=paste0(outputdir,"case1_spec_tad.bed"),sep='\t',row.names = F,col.names=F,quote=F)
write.table(tad_2_spec_10k,file=paste0(outputdir,"case2_spec_tad.bed"),sep='\t',row.names = F,col.names=F,quote=F)