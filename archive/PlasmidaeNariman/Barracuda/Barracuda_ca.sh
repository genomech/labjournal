#!/bin/bash


Threads=10;
ErrorRate=0.1;
Seq1R="";
Seq1F="";
OutputR1="";
OutputR2="";
InputR1="";
InputR2="";

cutadapt -j $Threads -e $ErrorRate -g $Seq1R -g $Seq1F -G $Seq1R -G $Seq1F -O 15 -o $OutputR1 -p $OutputR2 $InputR1 $InputR2
