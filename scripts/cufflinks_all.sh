 #!/bin/bash

THREADS=12
SORTED_PATH="/dev/datasets/FairWind/_results/Andre/bam_sorted"
CUFFLINKS_PATH="/dev/datasets/FairWind/_results/Andre/cufflinks"

# Cufflinks

mkdir -p $CUFFLINKS_PATH

for var in '001_S1' '002_S2' '003_S3' '004_S4' '005_S5' '006_S6' '007_S7' '008_S8' '009_S9' '010_S10' '011_S11' '012_S12' 
do
mkdir -p $CUFFLINKS_PATH/$var;
cufflinks -p $THREADS --library-type fr-firststrand \
	-o $CUFFLINKS_PATH/$var \
	$SORTED_PATH/MB_AFR_"$var"_sorted.bam;
echo $var is cuffled.
done

md5deep -lr $CUFFLINKS_PATH/* > $CUFFLINKS_PATH/all.md5;
chmod 555 -R $CUFFLINKS_PATH;

echo Cufflinks are sealed.
