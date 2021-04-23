#!/bin/bash exec 3>&1 4>&2 trap 'exec 2>&4 1>&3' 0 1 2 3 exec 1>log.out 2>&1
# Settings
FASTQFOLDER=fastq_pass;
MEMLIMIT=50g; # Adjust this based on your hardware (Do not use the maximum available)
THREADS=50; # Adjust this based on your hardware (You can get number of available threads using "nproc" and might want to save one for running other stuff)
MINQUALITY=80; # Remove reads with a quality of less than X %
DATABASE_FASTA=/srv/MA/users/mha/databases/midas36/midas36_notax.fa; # Path to fasta database for mapping
DATABASE_TAX=/srv/MA/users/mha/databases/midas36/tax_complete_qiime.txt # Path to taxonomy of database set appropriate error handling
######
set -o errexit -o pipefail -o noclobber -o nounset
echoWithHeader() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}
#check if a folder is present and empty
checkFolder() {
  if [ -d $1 ]
      then
          echoWithHeader "A directory named '$1' already exists. Do you want to continue (y/n)?"
          read ANSWER
          if [ $ANSWER != "y" ]
              then
                echoWithHeader "Exiting script..."
                echo ""
                exit 0
          fi
      else
        mkdir -p $1
  fi
}

# Start workflow
DATABASE_NAME=$(echo "$DATABASE_FASTA" | grep -o "[^/]*$" | grep -o "^[^\.]*")
echo ""
echoWithHeader "Running Nanopore raw mapping script for live demultiplexed files (max threads $THREADS)"
checkFolder temp
#checkFolder output
####
mkdir -p temp/mapped/
mkdir -p temp/filtered/

# Cut barcodes and adapters off
#module load Porechop/0.2.3-foss-2018a-Python-3.6.4
#printf "\n#Demultiplexing into '/barcodes/' folder----------------------------\n\n"
#mkdir -p barcodes
#for dirpath in $FASTQFOLDER/barcode*/; do
#    barcodename=$(basename "$dirpath")
#    echo "Processing $barcodename"
#    porechop -i $dirpath -b barcodes/ --threads 40 --no_split --check_reads 1000
#done
listoffastq=$(find $FASTQFOLDER/barcode*/ -name '*.fastq' -printf '%h\n' | sort -u)
for dirpath in $listoffastq/; do
barcodename=$(basename "$dirpath")
if [ -s $FASTQFOLDER/concatenated/$barcodename.fastq ]; then echo "$barcodename has already been concatenated";  
else
    (
    mkdir -p $FASTQFOLDER/concatenated/
	  echo "Concatenating $barcodename"
      cat $dirpath/*.fastq | \
      sed 's/\(runid.*start_\)//' | \
      tr -d '[:blank:]' > $FASTQFOLDER/concatenated/$barcodename.fastq
    )
fi
done  

# Filter reads based on quality
#module load Filtlong/0.2.0-foss-2018a
#for f in $FASTQFOLDER/concatenated/*.fastq; do
#NAME=$(basename "$f" .fastq)
#if [ -s temp/filtered/$NAME.filtered.fastq ]; then echo "temp/filtered/$NAME.filtered.fastq has already been generated";  
#else
#filtlong --min_mean_q $MINQUALITY $f > \
#  temp/filtered/$NAME.filtered.fastq
#fi


for f in $FASTQFOLDER/concatenated/*.fastq; do
#for f in temp/filtered/*.fastq; do
NAME=$(basename "$f" .fastq)
if [ -s temp/$NAME.idmapped.txt ]; then echo "temp/$NAME.idmapped.txt has already been generated";  
else
module load Minimap2/2.15-foss-2018a
# Map reads to reference database
minimap2 -ax map-ont \
  --cap-sw-mem=$MEMLIMIT \
  -t $THREADS \
  --secondary=no \
  $DATABASE_FASTA \
  -K20M $f > \
  temp/mapped/$NAME.sam
module load SAMtools/1.10-foss-2018a
samtools view -F 256 \
  -F 4 \
  -F 2048 temp/mapped/$NAME.sam \
  -o temp/mapped/$NAME"_nodupes.sam"
# Create mapping overview
echoWithHeader "Creating mapping overview"
sed '/^@/ d' temp/mapped/$NAME"_nodupes.sam" | \
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
   sed 's/time=/\t/' | sed 's/Zflow_cell.*barcode=/\t/' > temp/$NAME.idmapped.txt
fi
done

echoWithHeader "Generating OTU/mapping table for loading into ampvis2"

# Load R as module
module load R/3.5.0-foss-2018a-X11-20180131

R --slave --args "$THREADS" "$DATABASE_TAX" "$DATABASE_NAME" << 'makeOTUtable'
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
               header = FALSE,
               sep = "\t",
               colClasses = "character",
               col.names = c("OTU", "tax"))

#read mappings
mappings<-data.frame(OTU=NULL,SeqID=NULL)
files<-list.files(path = "temp/", pattern = ".idmapped.txt")
for (file in files){
  mapping <- fread(
    paste0("temp/",file),
    header = FALSE,
    #select = c(3),
    #colClasses = "character",
	col.names = c("readID", "read_time", "add_info")
  ) %>% separate(add_info, c("barcode", "SAMflag", "OTU", "Qlen", "alnlen", "MapID", "NMtag", "alnscore", "MinimapID"), " ") %>% mutate(SeqID=sub(".idmapped.txt","",file))
  mappings<-rbind(mappings,mapping)
}

# Subset mappings based on 
mappings_s <- mappings %>% mutate(Qr = as.numeric(Qlen)/as.numeric(alnlen)) %>% 
subset(., Qr < 1.15 & Qr > 0.85)
mappings_to_otu <- mappings_s %>% 
select(c("SeqID", "OTU"))

# Write mappings out for further analysis in R
question <- askYesNo("Do you want to write out a detailed mapping file?", default=TRUE, prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
if (question == TRUE) {
  fwrite(mappings_s, "temp/mappings_filt15.txt", quote = F, sep = "\t", row.names = F, col.names = T)
}else if (question == FALSE) {
  print("Okay, we'll resume..")
}

# Define function for 
#join taxonomy with mapping
joined <- taxDB[mappings_to_otu, on = "OTU"]

#transform into "OTU table", where rows are OTU's, columns are sample AND taxonomy
BIOMotutable <- dcast(joined, OTU + tax ~ SeqID, fun.aggregate = length) %>% setDT()
BIOMotutable[,taxonomy := gsub("[a-z]__", "", tax)]
BIOMotutable[,tax := NULL]
BIOMotutable[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") := tstrsplit(taxonomy, ";", fixed=FALSE)]
BIOMotutable[,taxonomy := NULL]
#write out
fwrite(BIOMotutable, paste0("temp/otutable_", args[[3]], ".txt"), sep = "\t", col.names = TRUE, na = "NA", quote = FALSE)

makeOTUtable

duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithHeader "Done in: $duration, enjoy!"
