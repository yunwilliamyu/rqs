Copyright (c) 2013 Y. William Yu. Released under CC0 1.0 Universal.

This package naively implements the Read-Quality-Sparsifier described in:
Y. William Yu, Deniz Yorukoglu, and Bonnie Berger. "Traversing the k-mer
landscape of NGS read datasets for quality score sparsification." Research in
Computational Molecular Biology, p385-399, 2014.
http://rqs.csail.mit.edu/

Note that the code described/implemented in this package is not scalable to
full genome / large datasets due to high memory requirements. This package was
written as a proof-of-principle for the RECOMB proceedings above and will *NOT*
be actively maintained.

-----------------------------
Package contents:

	fastq2fsam.sh:  converts FASTQ files to a SAM-like format
	fsam2fastq.sh:  converts SAM-like files to FASTQ

Following programs take input and output in SAM-like format:

	generate_dict:  builds a dictionary of common k-mers from a corpus
		 sparsify:  uses dictionary to smooth quality values for high confidence
					calls as measured by k-mer Hamming distance.
		 theshold:  reduces all quality values to some cutoff Q

-----------------------------
Dependencies:

	GCC 4.7 (or another C++11 compliant compiler)
	BOOST Multiindex

Quickstart:

	make
	./generate_dict MINCOUNT dict.db *.sam
	./sparsify dict.db *.sam
	./threshold 'Q' *.sam.filtered

We also provide a testsuite/ directory with example FASTQ files to
demonstrate operation. To run the commented example script:

	cd testsuite/
	./run.sh
-----------------------------
The three main programs are as follows.

generate_dict:
Usage: ./generate_dict MINCOUNT output_file input_file(s)

	Counts the number of times 32-mers appear, and outputs it in a two column
	format, with the first column specifying the 32-mer and the second column
	the number of times it appears in the corpus.

sparsify:
Discards likely non-SNP quality scores for known reads.
Usage: ./sparsify database_file input_file(s)

	Input is assumed to be SAM file without headers. Because this program
	only modifies column 11 if column 10 contains a valid read, however,
	it may or may not work on SAM files with headers. This behaviour
	should *not* be depended on.

	Output will be modified SAM files named [input_file].filtered,
	identical to the original except in column 11.

	The quality of nearly all bases corresponding to a 32-mer listed in the 1st
	column of [database_file] will be set to '~'. Correspondance shall be
	defined as a Hamming distance of less than or equal to 1 in the 32-mer.
	Note that a 32-mer can correspond to multiple 32-mers in the database.
	The exceptions to the setting of quality values above will be all bases
	that are different to any of the corresponding 32-mers in the database,
	which will retain their original quality value data.

	Example w/ 8-mers: suppose "AAAAAAAA" and "TAAAAAAT" are in the database,
	and we have the 8-mer "TAAAAAAA" in the input_file with original qualities
	"ABCDEFGH". Then the new quality values will be "A~~~~~~H".

	Note that for reads longer than 32 bases, filter will only check 32-mers
	for matches at most once every 16 bases, and once more at the end.

threshold:
Changes the maximum quality value of '~' to 'Q'
Usage: ./threshold 'Q' input_file(s)

	Will output a modified SAM file named '[input_file].reduced' with all
	instances of '~' in column 11 replaced with with character specified in
	[Q].

	In the previous example, the new quality values from the filtered SAM file
	would then be "AIIIIIIH" if Q=I. 
-----------------------------

For more resources and larger test cases, see http://rqs.csail.mit.edu/
