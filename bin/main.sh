#!/bin/bash
source activate motif 

## 1. get fasta
for i in *10.bed;do echo $i ;getFasta.sh <(cat $i| grep -E "chr([0-9,X,Y]+)") hg19 > ${i/bed/fa} & sleep 1; done

## 2. call FIMO/AME
for i in *10.fa; do echo $i; callAME.sh $i & sleep 1; done 

## 3. parse/gather AME's results
parseAME.sh

## 4. add tfclass annotations


############################################################
## use all di-nuc shuffle as background
############################################################






