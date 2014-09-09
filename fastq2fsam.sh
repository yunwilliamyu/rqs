#!/bin/bash
if [ $# -eq 0 ]
then
	echo Converts FASTQ files to a SAM-like format.
	echo Usage: $0 file1.fastq file2.fastq ...
	echo -e '\t'Generates file1.fastq.fsam file2.fastq.fsam ...
	echo -e '\t'Output files mimic the tab separated format of SAM, with read name is column 1, read sequence in column 10, and read quality in column 11.
	echo -e '\t'This is used to process files for input into the RQS package as generate_dict, sparsify, and threshold assume SAM-like input.
fi

for file in $@
do
	echo $file '-->' $file.fsam
	cat $file | sed '$!N;s/\n/\t/' | sed '$!N;s/\n/\t/' | sed 's/^@//' | \
		awk 'BEGIN { OFS=FS="\t" }; {print $1, "+", "+", "+", "+", "+", "+", "+", "+", $2, $4, "+"}' > \
		$file.fsam
done
