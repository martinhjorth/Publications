#!/bin/bash
module load Minimap2/2.17-foss-2018a
#MINIMAPDIR=map_db_test/minimap2/minimap2;

## Reads
#SIMREADS=map_db_test/simulated/badread_clost/guppy235_called/ESV13758_badread_guppy235.fastq;
#SIMREADS=seqdata/nanopore/badread/2019-10-10/v5/*.fastq;
#SIMREADS=/srv/MA/users/mha/Onsite_paper/seqdata/nanopore/badread/clostridiales-2k-906-meanID.fastq;
SIMREADS=/srv/MA/users/mha/Onsite_paper/seqdata/nanopore/2020-09-29_clostridiales-badread/fastq/clostridiales_badread_sim_v5.fastq;
#REALREADS=map_db_test/clostridiales_esv.fa;
REALREADS=/srv/MA/users/mha/Onsite_paper/map_db_test/clostridiales_esv_midas35.fa;

## Databases
#REF=db/midas33/midas33_single.fa;
#SPECIES=db/midas33-mod/blautia-midas33-species-filt.fa;
#GENUS=db/midas33-mod/blautia-midas33-genus-filt.fa;
#FAMILY=db/midas33-mod/blautia-midas33-family-filt.fa;
#ORDER=db/midas33-mod/blautia-midas33-order-filt.fa;
#CLASS=db/midas33-mod/blautia-midas33-class-filt.fa;
#PHYLUM=db/midas33-mod/blautia-midas33-phylum-filt.fa;
#DOMAIN=db/midas33-mod/blautia-midas33-domain-filt.fa;

REF=/srv/MA/users/mha/databases/midas37/midas37_single.fa;
SPECIES=/srv/MA/users/mha/Onsite_paper/db/midas37-mod/blautia-midas37-species-filt.fa
GENUS=/srv/MA/users/mha/Onsite_paper/db/midas37-mod/blautia-midas37-genus-filt.fa;
FAMILY=/srv/MA/users/mha/Onsite_paper/db/midas37-mod/blautia-midas37-family-filt.fa;
ORDER=/srv/MA/users/mha/Onsite_paper/db/midas37-mod/blautia-midas37-order-filt.fa;
CLASS=/srv/MA/users/mha/Onsite_paper/db/midas37-mod/blautia-midas37-class-filt.fa;
PHYLUM=/srv/MA/users/mha/Onsite_paper/db/midas37-mod/blautia-midas37-phylum-filt.fa;
DOMAIN=/srv/MA/users/mha/Onsite_paper/db/midas37-mod/blautia-midas37-domain-filt.fa;

## Function for processing SAM files
processSAM() {
awk '!/^@/' | awk -F "\t" '$2 == 0{
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
  }'
}


### Reference ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $REF $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_ref.txt
##### Real
minimap2 -ax sr --MD -t 20 $REF $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_ref.txt


### Species ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $SPECIES $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_filt_species.txt
##### Real
minimap2 -ax sr --MD -t 20 $SPECIES $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_filt_species.txt


### Genus ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $GENUS $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_filt_genus.txt
##### Real
minimap2 -ax sr  --MD -t 20 $GENUS $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_filt_genus.txt


### Family ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $FAMILY $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_filt_family.txt
##### Real
minimap2 -ax map-ont --MD -t 20 $FAMILY $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_filt_family.txt


### Order ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $ORDER $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_filt_order.txt
##### Real
minimap2 -ax map-ont --MD -t 20 $ORDER $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_filt_order.txt


### Class ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $CLASS $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_filt_class.txt
##### Real
minimap2 -ax map-ont --MD -t 20 $CLASS $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_filt_class.txt


### Phylum ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $PHYLUM $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_filt_phylum.txt
##### Real
minimap2 -ax map-ont --MD -t 20 $PHYLUM $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_filt_phylum.txt


### Domain ###
##### Simulated
minimap2 -ax map-ont --MD --secondary=no -t 20 $DOMAIN $SIMREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_badread_minimap_vs_midas37_filt_domain.txt
##### Real
minimap2 -ax map-ont --MD -t 20 $DOMAIN $REALREADS | sed '/^@/ d' | processSAM > map_db_test/mappings/clost_real_minimap_vs_midas37_filt_domain.txt