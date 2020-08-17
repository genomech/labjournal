#!/bin/bash

# Global
threads=10;
genome_fa="/dev/datasets/FairWind/_db/hg19/hg19.fa";
exome_bed="/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.sorted.bed";
min_coverage=10;
max_coverage=200;

# FastQC Check
fastqc -o $output_dir -t $threads $R1_fq;

# BWA Align
bwa mem -t $threads -v 1 $genome_fa $R1_fq $R2_fq | samtools view -bS -@ $threads -O BAM - > $output_bam;

# BAM Statistics
samtools stats -@ $threads --reference $genome_fa $file_bam > $stat_txt;

# IGV Prepare
PicardCommandLine SortSam SO=coordinate I=$input_bam O=$output_bam;
samtools index $output_bam;

# Duplicate Remove
PicardCommandLine SortSam SO=queryname I=$input_bam O=$output_bam_temp0;
PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true M=$picard_metrics_txt I=$output_bam_temp0 O=$dupless_bam;
python3 $daemonic_dir/strandless.py -f BAM -i $dupless_bam -o $strandless_bam -m $strandless_metrics_txt;
rm -f $output_bam_temp0;

# Filter Call
freebayes -0 --min-coverage $min_coverage --max-coverage $max_coverage -f $genome_fa -t $exome_bed -b $input_bam | vcflib vcfallelicprimitives > $output_vcf;

# Annotation
perl $annovar_dir/table_annovar.pl $input_vcf $annovar_dir/humandb -buildver $genome_assembly -protocol knownGene,ensGene,refGene,abraom,AFR.sites.2015_08,ALL.sites.2015_08,AMR.sites.2015_08,ASN.sites.2012_04,avgwas_20150121,avsift,avsnp150,cadd13,cg69,clinvar_20190305,cosmic70,dann,dbnsfp35c,dbscsnv11,EAS.sites.2015_08,eigen,esp6500_all,EUR.sites.2015_08,exac03,fathmm,gene4denovo201907,gerp++,gme,gnomad211_genome,gwava,hrcr1,icgc21,intervar_20180118,kaviar_20150923,ljb26_all,mcap13,mitimpact24,MT_ensGene,nci60,popfreq_all_20150413,regsnpintron,revel,SAS.sites.2015_08,snp142 --operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput --thread $threads;
python3 $daemonic_dir/annofit.py $input_vcf."$genome_assembly"_multianno.txt $output_csv;
