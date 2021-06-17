#!/bin/bash exec 3>&1 4>&2 trap 'exec 2>&4 1>&3' 0 1 2 3 exec 1>log.out 2>&1
# Settings
FASTQFOLDER=$1;
APPEND=$2;
TODAY=$(date '+%Y-%m-%d')
OUT_DIR=${TODAY}"_processed"
MEMLIMIT=50g; # Adjust this based on your hardware (Do not use the maximum available)
THREADS=50; # Adjust this based on your hardware (You can get number of available threads using "nproc" and might want to save one for running other stuff)
MINLENGTH=2000; # Remove reads with a quality of less than X %
DATABASE_FASTA=db/midas-genomedb/HQ_sp_581_concat-w-medium-tricho.fa; # Path to fasta database for mapping
DATABASE_TAX=db/midas-genomedb/midas-genomedb-tax_w-trichococcus.txt # Path to taxonomy of database set appropriate error handling
######
set -o errexit -o pipefail -o noclobber -o nounset
echoWithHeader() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}


# Start workflow
DATABASE_NAME="midas_gDB-w-CMStricho"
echo ""
echoWithHeader "Running Nanopore raw mapping script (max threads $THREADS)"

####
mkdir ${OUT_DIR}
mkdir -p ${OUT_DIR}/temp/mapped/
mkdir -p ${OUT_DIR}/temp/filtered/
###

for dirpath in $FASTQFOLDER/barcode*/; do
barcodename=$(basename "$dirpath")
if [ -s ${OUT_DIR}/concatenated/$barcodename.fastq ]; then echo "$barcodename has already been concatenated";  
else
    (
    mkdir -p ${OUT_DIR}/concatenated/
	count=`ls -1 $dirpath/*.fastq 2>/dev/null | wc -l`
	if [ $count != 0 ]; then
	  echo "Concatenating $barcodename"
      cat $dirpath/*.fastq | \
      sed 's/\(runid.*start_\)//' | \
      tr -d '[:blank:]' > ${OUT_DIR}/concatenated/$barcodename.fastq
	else
	  continue
	fi
    )
fi
done  

# Filter reads based on length
module load Filtlong/0.2.0-foss-2018a
for f in ${OUT_DIR}/concatenated/*.fastq; do
NAME=$(basename "$f" .fastq)
if [ -s ${OUT_DIR}/temp/filtered/$NAME.filtered.fastq ]; then echo "temp/filtered/$NAME.filtered.fastq has already been generated";  
else
filtlong --min_length $MINLENGTH $f > \
  ${OUT_DIR}/temp/filtered/$NAME.filtered.fastq
fi
done

for f in ${OUT_DIR}/temp/filtered/*.fastq; do

NAME=$(basename "$f" .fastq)
if [ -s ${OUT_DIR}/temp/$NAME.idmapped.txt ]; then echo "temp/$NAME.idmapped.txt has already been generated";  
else
module load Minimap2/2.17-foss-2018a

# Map reads to reference database
minimap2 -ax map-ont \
  --cap-sw-mem=$MEMLIMIT \
  -t $THREADS \
  --secondary=no \
  $DATABASE_FASTA \
  -K20M $f > \
  ${OUT_DIR}/temp/mapped/$NAME.sam
module load SAMtools/1.10-foss-2018a
samtools view -F 256 \
  -F 4 \
  -F 2048 ${OUT_DIR}/temp/mapped/$NAME.sam \
  -o ${OUT_DIR}/temp/mapped/$NAME"_nodupes.sam"
  
# Create mapping overview
echoWithHeader "Creating mapping overview"
sed '/^@/ d' ${OUT_DIR}/temp/mapped/$NAME"_nodupes.sam" | \
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
  }' | \
   sed 's/time=/\t/' | sed 's/Zbarcode=/\t/' > ${OUT_DIR}/temp/$NAME.idmapped.txt
fi
done

echoWithHeader "Generating OTU/mapping table for loading into ampvis2"

# Load R as module
module load R/3.5.0-foss-2018a-X11-20180131

R --slave --args "$THREADS" "$DATABASE_TAX" "$DATABASE_NAME" "$OUT_DIR" << 'makeOTUtable'
#extract passed args from shell script
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({
  if(!require("data.table")) {
    install.packages("data.table")
    require("data.table")
  }
  if(!require("dplyr")) {
    install.packages("dplyr")
    require("dplyr")
  }
  if(!require("tidyr")) {
    install.packages("tidyr")
    require("tidyr")
  }
})

setDTthreads(as.integer(args[[1]]))

#read taxonomy
taxDB <- fread(args[[2]],
               header = TRUE,
               sep = "\t",
               colClasses = "character")

#read mappings
mappings<-data.frame(OTU=NULL,SeqID=NULL)
files<-list.files(path = paste0(args[[4]],"/temp/"), pattern = ".idmapped.txt")
for (file in files){
  mapping <- fread(
    paste0(args[[4]],"/temp/",file),
    header = FALSE,
	col.names = c("readID", "read_time", "add_info")
  ) %>% separate(add_info, c("barcode", "SAMflag", "OTU", "Qlen", "alnlen", "MapID", "NMtag", "alnscore", "MinimapID"), " ") %>% mutate(SeqID=sub(".idmapped.txt","",file)) %>% separate(OTU, c("MAG", "contig"), "~")
  mappings<-rbind(mappings,mapping)
}

# Subset mappings based on
mappings_s <- mappings %>% mutate(Qr = as.numeric(Qlen)/as.numeric(alnlen))#Filtering step can be added here!
mappings_to_otu <- mappings_s %>%
select(c("SeqID", "MAG"))

# Write mappings out for further analysis in R
question <- askYesNo("Do you want to write out a detailed mapping file?", default=TRUE, prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
if (question == TRUE) {
  fwrite(mappings_s, paste0(args[[4]],"/",format(Sys.time(), "%Y-%m-%d"), "_mappings_metag_", args[[3]], "_2kbpfilt.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
}else if (question == FALSE) {
  print("Okay, we'll resume..")
}

# Define function for
#join taxonomy with mapping
joined <- taxDB[mappings_to_otu, on = "MAG"]

#transform into "OTU table", where rows are OTU's, columns are sample AND taxonomy
BIOMotutable <- dcast(joined, MAG + Kingdom + Phylum + Class + Order + Family + Genus + Species ~ SeqID, fun.aggregate = length) %>% setDT()

#write out
fwrite(BIOMotutable, paste0(args[[4]],"/", format(Sys.time(), "%Y-%m-%d"), "_otutable_metag-tricho_", args[[3]], "_2kbpfilt.txt"), sep = "\t", col.names = TRUE, na = "NA", quote = FALSE)

makeOTUtable

duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithHeader "Done in: $duration, enjoy!"
