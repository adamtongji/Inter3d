args=commandArgs(T)

ctcf_gain_annot<-args[1]
ctcf_loss_annot<-args[2]
up_gene<-args[3]
down_gene<-args[4]
outdir<-args[5]
peak_input<-args[6]
# ctcf_gain_annot<-"./output/tad_case1_ctcf_annot.txt"
# ctcf_loss_annot<-"./output/tad_case2_ctcf_annot.txt"

# up_gene<-"testdata/DEG/up_gene.txt"
# down_gene<-"testdata/DEG/down_gene.txt"
# outdir<-"./output/"

ctcf_gain <- read.delim(ctcf_gain_annot,sep='\t',header=T)
ctcf_gain_gene<-unique(subset(ctcf_gain, abs(ctcf_gain$distanceToTSS)<100000)$SYMBOL)

ctcf_loss <- read.delim(ctcf_loss_annot,sep='\t',header=T)
ctcf_loss_gene<-unique(subset(ctcf_loss, abs(ctcf_loss$distanceToTSS)<100000)$SYMBOL)

# DEG must GENE OFFICIAL SYMBOL
rna_seq_up<-read.table(up_gene)
rna_seq_down<- read.table(down_gene)

# ===========================================================
# STEP 2: Identify candidate genes regulated by CTCF binding.
# ===========================================================

# Mechanism 1

case1_spec_gene<-unique(rna_seq_up[,1]) # ctcf loss
case2_spec_gene<-unique(rna_seq_down[,1]) # ctcf gain

case1_spec_gene_withctcf<-intersect(case1_spec_gene,ctcf_loss_gene)
case2_spec_gene_withctcf<-intersect(case2_spec_gene,ctcf_gain_gene)

ctcf_mech1_case1_genespec<-subset(ctcf_loss,SYMBOL %in% case1_spec_gene_withctcf)[,c(1:4,13,15:18,21:24)]


ctcf_mech1_case2_genespec<-subset(ctcf_gain,SYMBOL %in% case2_spec_gene_withctcf)[,c(1:4,6,8:11,14:17)]
#rnaseq_tab<-rbind(rna_seq_up[,c(3:6,10,14)],rna_seq_down[,c(3:6,10,14)])
#rnaseq_tab$SYMBOL=c(as.vector(rna_seq_up[,1]),as.vector(rna_seq_down[,1]))

# ctcf_case2_genespec_tab<-merge(ctcf_case2_genespec,rnaseq_tab,by="SYMBOL")
# ctcf_case1_genespec_tab<-merge(ctcf_case1_genespec,rnaseq_tab,by="SYMBOL")

write.table(ctcf_mech1_case2_genespec,file = paste0(outdir,"/ctcf_gain_case1_low_gene.xls"),sep='\t',quote=F,row.names = F,col.names = T)

write.table(ctcf_mech1_case1_genespec,file = paste0(outdir,"/ctcf_loss_case1_high_gene.xls"),sep='\t',quote=F,row.names = F,col.names = T)

# Mechanism 2

case1_regu_gene_withctcf<-intersect(case1_spec_gene,ctcf_gain_gene)
case2_regu_gene_withctcf<-intersect(case2_spec_gene,ctcf_loss_gene)

ctcf_mech2_case1_genespec<-subset(ctcf_gain,SYMBOL %in% case1_regu_gene_withctcf)[,c(1:4,6,8:11,14:17)]

ctcf_mech2_case2_genespec<-subset(ctcf_loss,SYMBOL %in% case2_regu_gene_withctcf)[,c(1:4,13,15:18,21:24)]

# rnaseq_tab<-rbind(rna_seq_up[,c(3:6,10,14)],rna_seq_down[,c(3:6,10,14)])
# rnaseq_tab$SYMBOL=c(as.vector(rna_seq_up[,1]),as.vector(rna_seq_down[,1]))

# ctcf_case2_generegu_tab<-merge(ctcf_case2_generegu,rnaseq_tab,by="SYMBOL")
# ctcf_case1_generegu_tab<-merge(ctcf_case1_generegu,rnaseq_tab,by="SYMBOL")

write.table(ctcf_mech2_case1_genespec,file = paste0(outdir,"/ctcf_gain_case1_high_gene.xls"),sep='\t',quote=F,row.names = F,col.names = T)
write.table(ctcf_mech2_case2_genespec,file = paste0(outdir,"/ctcf_loss_case1_low_gene.xls"),sep='\t',quote=F,row.names = F,col.names = T)



# ===========================================================
# STEP 3: Identify candidate genes regulated by CTCF binding.
# ===========================================================
share_peaks<-read.delim(peak_input,sep = "\t",header=T)

mech1_case1_genes<-as.vector(ctcf_mech1_case1_genespec$SYMBOL)

select_mech1_case1_cre<-c()
for (gene in mech1_case1_genes){
  relate_peaks<-subset(share_peaks, SYMBOL==gene)
  if (nrow(relate_peaks)>0){
    dist1 = subset(ctcf_mech1_case1_genespec,SYMBOL==gene)$distanceToTSS
    for (line in 1:nrow(relate_peaks)){
      dist2 = relate_peaks[line,"distanceToTSS"]
      if (abs(dist2)>abs(dist1) &  as.numeric(dist2)*as.numeric(dist1)>0){
        select_mech1_case1_cre<-rbind(select_mech1_case1_cre,relate_peaks[line,])
      }
    }
  }
}
select_mech1_case1_cre<-select_mech1_case1_cre[!duplicated(select_mech1_case1_cre[,2]),]
write.table(select_mech1_case1_cre,paste0(outdir,"/ctcf_gain_case1_low_gene_CRE.xls"),sep='\t',row.names = F,col.names = T,quote=F)

mech1_case2_genes<-as.vector(ctcf_mech1_case2_genespec$SYMBOL)

select_mech1_case2_cre<-c()
for (gene in mech1_case2_genes){
  relate_peaks<-subset(share_peaks, SYMBOL==gene)
  if (nrow(relate_peaks)>0){
    dist1 = subset(ctcf_mech1_case2_genespec,SYMBOL==gene)$distanceToTSS
    for (line in 1:nrow(relate_peaks)){
      dist2 = relate_peaks[line,"distanceToTSS"]
      if (abs(dist2)>abs(dist1) &  as.numeric(dist2)*as.numeric(dist1)>0){
        select_mech1_case2_cre<-rbind(select_mech1_case2_cre,relate_peaks[line,])
      }
    }
  }
}
select_mech1_case2_cre<-select_mech1_case2_cre[!duplicated(select_mech1_case2_cre[,2]),]
write.table(select_mech1_case2_cre,paste0(outdir,"/ctcf_loss_case1_high_gene_CRE.xls"),sep='\t',row.names = F,col.names = T,quote=F)

mech2_case1_genes<-as.vector(ctcf_mech2_case1_genespec$SYMBOL)

select_mech2_case1_cre<-c()
for (gene in mech2_case1_genes){
  relate_peaks<-subset(share_peaks, SYMBOL==gene)
  if (nrow(relate_peaks)>0){
    dist1 = subset(ctcf_mech2_case1_genespec,SYMBOL==gene)$distanceToTSS
    for (line in 1:nrow(relate_peaks)){
      dist2 = relate_peaks[line,"distanceToTSS"]
      if (abs(dist2)>abs(dist1) &  as.numeric(dist2)*as.numeric(dist1)>0){
        select_mech2_case1_cre<-rbind(select_mech2_case1_cre,relate_peaks[line,])
      }
    }
  }
}
write.table(select_mech2_case1_cre,paste0(outdir,"/ctcf_gain_case1_high_gene_CRE.xls"),sep='\t',row.names = F,col.names = T,quote=F)
select_mech2_case1_cre<-select_mech2_case1_cre[!duplicated(select_mech2_case1_cre[,2]),]
mech2_case2_genes<-as.vector(ctcf_mech2_case2_genespec$SYMBOL)

select_mech2_case2_cre<-c()
for (gene in mech2_case2_genes){
  relate_peaks<-subset(share_peaks, SYMBOL==gene)
  if (nrow(relate_peaks)>0){
    dist1 = subset(ctcf_mech2_case2_genespec,SYMBOL==gene)$distanceToTSS
    for (line in 1:nrow(relate_peaks)){
      dist2 = relate_peaks[line,"distanceToTSS"]
      if (abs(dist2)>abs(dist1) &  as.numeric(dist2)*as.numeric(dist1)>0){
        select_mech2_case2_cre<-rbind(select_mech2_case2_cre,relate_peaks[line,])
      }
    }
  }
}
select_mech2_case2_cre<-select_mech2_case2_cre[!duplicated(select_mech2_case2_cre[,2]),]
write.table(select_mech2_case2_cre,paste0(outdir,"/ctcf_loss_case1_low_gene_CRE.xls"),sep='\t',row.names = F,col.names = T,quote=F)



