###################################################################################################
#
#  Version 1.1
#
#  This workflow script uses cutadapt to trim adapter sequences from reads. Specifically it creates
#  two folders (trimmed and untrimmed), concatenates fastq files in dir barcodes, looking for porechop barcode names "BC01-12",
#  trims for the V1 and V3 primer sequences with an error rate of up to 15% (Nanopore read error),
#  outputs trimmed reads to the 'trimmed folder' and untrimmed to the 'untrimmed' folder. The
#  reverse complement sequences are then searched for in the untrimmed files and output to the
#  trimmed folder.
#  Run the script from the 'pass' folder to access all the 'barcodeXX' folders.
#
#  Martin Andersen, 2018
#
###################################################################################################

clear
echo ""
echo "Running cutadapt workflow version 1.1"
date

mkdir trimmed
mkdir untrimmed

for i in $(seq -w 01 12);
do
	cat barcode"$i"/*.fastq > barcode"$i"/barcode"$i"_combined.fastq
	cutadapt -g AGAGTTTGATCCTGGCTCAG...CCAGCAGCCGCGGTAAT -e 0.15 --untrimmed-output untrimmed/barcode"$i"_untrimmed.fastq -o trimmed/barcode"$i"_trimmed1.fastq barcode"$i"/barcode"$i"_combined.fastq
	cutadapt -g ATTACCGCGGCTGCTGG...CTGAGCCAGGATCAAACTCT -e 0.15 --discard-untrimmed -o trimmed/barcode"$i"_trimmed2.fastq untrimmed/barcode"$i"_untrimmed.fastq
	#sed -i "s/sampleid=LIB-MHA-62-A-2/barcode=barcode$i/" trimmed/barcode"$i"_trimmed*.fastq
done

echo ""
echo "Concatenating trimmed reads into 'np_all_combined.fastq' file"
cat trimmed/*.fastq > np_all_gup315_trimmed.fastq

notrimmed=$(grep '^@.*' np_all_gup315_trimmed.fastq | wc -l | xargs)
original=$(find barcode"$i"/barcode"$i"_combined.fastq -name -exec cat {} \; | grep '^@.*' | wc -l | xargs)

percentage=$(echo "(($notrimmed / $original) * 100)" | bc -l)
echo "Done! You have $notrimmed trimmed reads"
printf "Or %.2f %% of the barcoded reads\n" $percentage
