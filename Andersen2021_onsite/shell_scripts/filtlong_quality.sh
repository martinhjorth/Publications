readinput=np_all_gup315_trimmed.fastq;
quality=80;
flversion=$(filtlong --version)

echo "Running $flversion"
echo ""
filtlong --min_mean_q $quality $readinput > np_all_gup315_trim_filt.fastq

inputno=$(grep '^@.*' $readinput | wc -l | xargs)
outputno=$(grep '^@.*' np_all_gup315_trim_filt.fastq | wc -l | xargs)

percentage=$(echo "(($outputno / $inputno) * 100)" | bc -l)
echo "Done! You have $outputno trimmed reads"
printf "Or %.2f %% of the input reads\n" $percentage
