 #!/bin/bash

THREADS=12
BAM_PATH="/dev/datasets/FairWind/_results/Fatima/bam"
SORTED_PATH="/dev/datasets/FairWind/_results/Fatima/sorted"
CUFFLINKS_PATH="/dev/datasets/FairWind/_results/Fatima/cufflinks"

# Досортировка

#mkdir -p $SORTED_PATH

#for var in '001_S13' '002_S14' '003_S15' '004_S16' '005_S17' '006_S18' '007_S19' '008_S20' '009_S21' '010_S22'
#do
#samtools sort -@ $THREADS -O BAM $BAM_PATH/MB_FQ_"$var".bam > $SORTED_PATH/MB_FQ_"$var"_sorted.bam;
#echo $var is sorted.
#done

#md5sum $SORTED_PATH/*.bam > $SORTED_PATH/all.md5;
#chmod 555 -R $SORTED_PATH;

#echo Sorted are sealed.

# Cufflinks

mkdir -p $CUFFLINKS_PATH

for var in '006_S18' '010_S22' 
do
mkdir -p $CUFFLINKS_PATH/$var;
cufflinks -p $THREADS --library-type fr-firststrand \
	-o $CUFFLINKS_PATH/$var \
	$SORTED_PATH/MB_FQ_"$var"_sorted.bam;
echo $var is cuffled.
done

md5deep -lr $CUFFLINKS_PATH/* > $CUFFLINKS_PATH/all.md5;
chmod 555 -R $CUFFLINKS_PATH;

echo Cufflinks are sealed.
