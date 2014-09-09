#!/bin/bash
for file in example?.fastq
do
	echo Cleaning up files generated from $file
	rm -rf	$file.fsam \
			$file.fsam.filtered \
			$file.fsam.filtered.reduced \
			$file.fsam.filtered.reduced.qual \
			$file.fsam.filtered.reduced.qual.bz2 \
			$file.fsam.filtered.reduced.fastq
done
echo "Removing generated dictionary dict.db"
rm -rf dict.db
