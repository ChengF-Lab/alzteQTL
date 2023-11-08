#!/usr/bin/bash

input="squire_map"
sample=`cat $1`

for sam in $sample
	do
		echo "START at date: "`date`
		
		TEcount --sortByPos --format BAM --mode multi -b $input/$sam\.bam --GTF hg38.refGene.gtf --TE GRCh38_GENCODE_rmsk_TE.gtf --project $sam -i 100

		echo "END at date: "`date`
done


