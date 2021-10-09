# step 1: Identify CTCF-mediated TAD alterations.

Rscript ../bin/tad_spec.r ./testdata/tad_boundary/TAD_genome1.bed ./testdata/tad_boundary/TAD_genome2.bed ./output/

intersectBed -a ./testdata/CTCF/ctcf1.bed.gz -b testdata/CTCF/ctcf2.bed.gz -v  >./output/ctcf1_spec.bed
intersectBed -a ./testdata/CTCF/ctcf2.bed.gz -b testdata/CTCF/ctcf1.bed.gz -v  >./output/ctcf2_spec.bed

intersectBed -a ./output/ctcf1_spec.bed -b ./output/case1_spec_tad.bed -wa > ./output/tad_case1_ctcf.bed
intersectBed -a ./output/ctcf2_spec.bed -b ./output/case2_spec_tad.bed -wa > ./output/tad_case2_ctcf.bed


echo "step 1 complete!"

# step 2: Identify candidate genes regulated by CTCF binding.

Rscript ../bin/annot_ctcf.r ./output/tad_case1_ctcf ./output/tad_case2_ctcf

