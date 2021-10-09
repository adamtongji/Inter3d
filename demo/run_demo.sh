# step 1: Identify CTCF-mediated TAD alterations.

Rscript ../bin/tad_spec.r ./testdata/tad_boundary/TAD_genome1.bed ./testdata/tad_boundary/TAD_genome2.bed ./output/

intersectBed -a ./testdata/CTCF/ctcf1.bed.gz -b testdata/CTCF/ctcf2.bed.gz -v  >./output/ctcf1_spec.bed
intersectBed -a ./testdata/CTCF/ctcf2.bed.gz -b testdata/CTCF/ctcf1.bed.gz -v  >./output/ctcf2_spec.bed

intersectBed -a ./output/ctcf1_spec.bed -b ./output/case1_spec_tad.bed -wa > ./output/tad_case1_ctcf.bed
intersectBed -a ./output/ctcf2_spec.bed -b ./output/case2_spec_tad.bed -wa > ./output/tad_case2_ctcf.bed


echo "step 1 complete!"

# step 2+3: Identify candidate genes regulated by CTCF binding. And then find Enhancer-Promoter pairs regulated by CTCF binding.

intersectBed -a testdata/atac_peaks/atac1.narrowPeak.gz -b testdata/atac_peaks/atac2.narrowPeak.gz -wa > ./output/share_peaks.bed

Rscript ../bin/annot_hg19_region.r ./output/tad_case1_ctcf ./output/tad_case2_ctcf ./output/share_peaks


Rscript ../bin/find_mech_genes_peak.r ./output/tad_case1_ctcf_annot.txt ./output/tad_case2_ctcf_annot.txt testdata/DEG/up_gene.txt testdata/DEG/down_gene.txt ./output/ ./output/share_peaks_annot.txt



