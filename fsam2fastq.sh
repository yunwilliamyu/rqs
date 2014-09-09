#!/bin/bash
if [ $# -eq 0 ]
then
	echo Converts SAM-like files to FASTQ format.
	echo Usage: $0 file1.fsam file2.fsam ...
	echo -e '\t'Generates file1.fsam.fastq file2.fsam.fastq ...
	echo -e '\t'Input files mimic the tab separated format of SAM, with read name is column 1, read sequence in column 10, and read quality in column 11.
	echo -e '\t'This is used to convert files converted by fastq2fsam back to FASTQ.
fi

for file in $@
do
	echo $file '-->' $file.fastq
	cat $1 | awk '{ print "@" $1 "\n" $10 "\n+\n" $11  }' > \
		$file.fastq
done
