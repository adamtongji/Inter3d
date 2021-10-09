ctcf1=$1
ctcf2=$2
tad1=$3
tad2=$4
atacpeak1=$5
atacpeak2=$6

outdir=$7
tmpdir=${outdir}/tmp

mkdir -p $tmpdir
intersectBed -a $ctcf1 -b $ctcf2 -v -f 0.1 > $tmpdir/ctcf1.bed
intersectBed -a $ctcf2 -b $ctcf1 -v -f 0.1 > $tmpdir/ctcf2.bed
