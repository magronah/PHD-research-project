#!/bin/bash
for x in *.fastq.gz
do 
	cutadapt -g ACTCCTACGGGAGGCAGCAG $x -o ${x%.fastq.gz}_trim.fastq.gz -q 30 -m 100
done
#echo

# echo mv *2.fastq.gz ~/Documents/data_PRJEB29421

