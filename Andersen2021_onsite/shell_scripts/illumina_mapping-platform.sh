READS=seqdata/illumina/il_midasSTD.fastq;
#DB=map_db_test/midas3/singleline_ESV.fa;
DB=db/midas31/MiDAS-31-single.fa;
mapper_path=/space/users/mha/Desktop/Onsite_paper/map_db_test/minimap2/minimap2;

$mapper_path -ax sr -t 20 --secondary=no --MD $DB $READS > map_db_test/mappings/illumina_midas31.sam

echo ""
echo "Mapped reads. Removing supplementary mappings based on SAM bit flags"
samtools view -F 256 -F 4 -F 2048 map_db_test/mappings/illumina_midas31.sam -o map_db_test/mappings/il_midas31_nodupes.sam
date

sed '/^@/ d' map_db_test/mappings/il_midas31_nodupes.sam | \
	awk '{
    for(i=1;i<=NF;i++){
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); mm = $i}
    }
    split($6, count, /[^0-9]+/);
    split($6, type, /[^A-Z]*/);
    for(i=1; i<= length(count)-1; i++){
      if(type[i + 1] ~ /[DIM]/){aln+=count[i]};
    }
    print $1, $2, $3, length($10), aln, (aln - mm)/aln, $12, $14, $20
    aln=0;
  }' > idmapped_ilv13.txt

awk -F"\t" -v OFS='\t' '{print $1, $3}' idmapped_ilv13.txt | sed 's/Zbarcode=/\t/' >  map_db_test/mappings/minimap_il_midas31_processed.txt
