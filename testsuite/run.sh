#!/bin/bash
echo First convert FASTQ files to a fake SAM format:
echo "../fastq2fsam.sh example.fastq example2.fastq"
../fastq2fsam.sh example?.fastq
echo

echo Generate our dictionary with read multiplicity r=2:
echo "../generate_dict 2 dict.db *.fsam"
../generate_dict 2 dict.db *.fsam
echo

echo In this example, we compress the corpus we used to generate the dictionary:
echo "../sparsify dict.db *.fsam"
../sparsify dict.db *.fsam
echo

echo After sparsification, we still need to cut-off high quality scores:
echo "Let's choose a threshold of 'I', or Phred quality 40 under most encodings"
echo "../threshold 'I' *.fsam.filtered"
../threshold 'I' *.fsam.filtered
echo

echo Compression is now done, but let\'s see how well we did using BZIP2:
echo "We'll cut out column 11 (the quality scores) and pipe it through BZIP2"
echo and then compute the bits per quality value needed.
for file in *.filtered.reduced
do
        cut -f11 $file > $file.qual
        orig_size=`wc -c < $file.qual`
        orig_lines=`wc -l < $file.qual`
        orig_size=`echo "$orig_size - $orig_lines" | bc`
        bzip2 -f $file.qual
        new_size=`wc -c < $file.qual.bz2`
        #rm $file.qual.bz2
        echo -e $file:'\t' `echo "scale=4; 1/( $orig_size / ( $new_size * 8)) " | bc` bits / quality score
done
echo

echo "Although we're basically done, you might want to convert your files back"
echo "from this fake SAM format to a FASTQ file:"
echo "../fsam2fastq.sh *.filtered.reduced"
../fsam2fastq.sh *.filtered.reduced
echo

echo "All done!"
echo
echo "If you want to clean up all the files generated in this example, just run"
echo "./cleanup.sh"
