#!/bin/bash
#######################################################
###
### input: http://tfclass.bioinf.med.uni-goettingen.de/suppl/tfclass.ttl.gz
### output genus.txt
### output subfamily.txt
### output family.txt
### output human_tf.txt
### output mouse_tf.txt 
###
### 
### Details: 
#####################################################

### genus - least level 
grep -E -A2 "^<http://sybig.de/tfclass#\d+.\d+.\d+.\d+.\d+>" tfclass.ttl |\
    grep -v "Genus" | \
    perl -pe 's/>(.*)\n/\t/mg' | \
    perl -ne '/(\d+.\d+.\d+.\d+.\d+)(.*)"(.*)"/ && print "$1\t$3\n"'  > genus.txt 

### subfamily -with seq 
grep -E -A3 "^<http://sybig.de/tfclass#\d+.\d+.\d+.\d+>" tfclass.ttl | \
    grep -v "Subfamily" | \
    perl -pe 's/>(.*)\n/>\t/mg' |\
    perl -pe 's/;(.*)\n/>\t/mg' |\
    perl -ne '/(\d+.\d+.\d+.\d+)>(.*)label "(.*)"(.*)sequence "(.*)"/ && print "$1\t$3\t$5\n"'  > subfamily.txt

### subfamily -all 
grep -E -A3 "^<http://sybig.de/tfclass#\d+.\d+.\d+.\d+>" tfclass.ttl | \
    grep -v "Subfamily" | \
    perl -pe 's/>(.*)\n/>\t/mg' |\
    perl -pe 's/;(.*)\n/>\t/mg' |\
    perl -ne '/(\d+.\d+.\d+.\d+)>(.*)label "(.*)"/ && print "$1\t$3\n"' |\
    perl -pe 's/"(.*)"/\t/mg'> subfamily_all.txt


### family

grep -E -A2 "^<http://sybig.de/tfclass#\d+.\d+.\d+>" tfclass.ttl | \
    grep -v "Family" | \
    perl -pe 's/>(.*)\n/\t/mg' |\
    perl -ne '/(\d+.\d+.\d+)(.*)"(.*)"/&&print "$1\t$3\n"' > family.txt

### for human : l1- get name +id l2- put into one line; l3 - extract name +id
perl -ne "/^<(.*)Homo_sapiens(.*)>|(:xref(.*)ENSG(\d+))/ && print " tfclass.ttl|\
    perl -pe 's/(Homo_sapiens_(.*))\n/$1\t/mg'|\
    perl -ne '/Homo_sapiens_(.*)(>(.*))(ENSG(\d+))/ &&  print "$1 \t $4 \n"' |\
    sort | uniq > human_tf.txt

# deal with irregular profile 
awk -v FS='\t' '(NF!=2)' ./human_tf.txt | perl -ne '/Homo_sapiens_(\S+)(\s+)(ENSG(\d+))/ && print "$1\t$3\n"' >tmp
cat <(awk -v FS='\t' '(NF==2)' ./human_tf.txt) tmp |sort | uniq  > human_tf.txt



### for mouse: l1- get name +id l2- put into one line; l3 - extract name +id 
perl -ne "/^<(.*)Mus_musculus(.*)>|(:xref(.*)ENSMUSG(\d+))/ && print " tfclass.ttl| \ 
    perl -pe 's/(Mus_musculus_(.*))\n/$1\t/mg'| \ 
    perl  -ne '/Mus_musculus_(.*)(>(.*))(ENSMUSG(\d+))/ &&  print "$1 \t $4 \n"'| \ 
    sort | uniq > mouse_tf.txt

cat <(awk -v FS='\t' '(NF!=2)' ./mouse_tf.txt | perl -ne '/Mus_musculus_(\S+)(\s+)(ENSMUSG(\d+))/ && print "$1\t$3\n"' )
    


### check TF numbers     
wc -l mouse_tf.txt # 1244
wc -l human_tf.txt # 1521 
