#!/bin/sh



###################################
#
# This script sits in the uparse directory and executes from there. 
# When 
#
###################################


# first, set up macqiime so we can run these commands in a subshell
# this only necessary on a mac. On ACISS, load QIIME and UPARSE instead. 
# source /macqiime/configs/bash_profile.txt


########################
# Set some directory shortcuts
#   !! This assumes a default macqiime 1.8 install and the same directory architecture. 

# export QIIME_DIR=/macqiime
# export reference_seqs=$QIIME_DIR/greengenes/gg_13_8_otus/rep_set/97_otus.fasta
# export reference_tax=$QIIME_DIR/greengenes/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
# export reference_tree=$QIIME_DIR/greengenes/gg_13_8_otus/trees/97_otus.tree
# export COMBINED=../rawData/combinedRuns
# export UNCOMBINED=../rawData



#############  commands to run interactive node on aciss cluster

qsub -q fatnodes -I

module load usearch 
module load fastx_toolkit/0.0.13
module load qiime/1.8.0  
printenv PATH # check to make sure lots of qiime dependencies are loaded. 

# Qiime has so many dependencies, it is difficult to load them all.
# So might have to reload module version. This makes no sense. 
source  /usr/local/packages/Modules/setmodule  3.2.9
source  /usr/local/packages/Modules/setmodule  3.2.10

################


##### Check quality - this goes into R script Pickle2014/R/quality/
# split_libraries_fastq.py -v -q 0 -i raw1/r1readTRIM.fastq -b raw1/barcodesRenamed.fastq -o splitLib1F/ -m map.txt --barcode_type 16 
# split_libraries_fastq.py -v -q 0 -i raw2/r1readTRIM.fastq -b raw2/barcodesRenamed.fastq -o splitLib2F/ -m map.txt --barcode_type 16 
# split_libraries_fastq.py -v -q 0 -i raw1/r2readTRIM.fastq -b raw1/barcodesRenamed.fastq -o splitLib1R/ -m map.txt --barcode_type 16 
# split_libraries_fastq.py -v -q 0 -i raw2/r2readTRIM.fastq -b raw2/barcodesRenamed.fastq -o splitLib2R/ -m map.txt --barcode_type 16 

# Ended up combining forward reads and just analyzing these instead of joining.  
cat raw1/r1readTRIM.fastq raw2/r1readTRIM.fastq > seqs.fastq
cat raw1/barcodesRenamed.fastq raw2/barcodesRenamed.fastq > barcodes.fastq
split_libraries_fastq.py -v -q 0 --store_demultiplexed_fastq -i seqs.fastq -b barcodes.fastq -o splitLib/ -m map.txt --barcode_type 16 # -n 300


# trim to 250 length. This is not straightforward or well documented in UPARSE, so 
# farm out to fastx. 
fastx_trimmer -l 250 -i splitLib/seqs.fastq -o splitLib/seqs.trimmed.fastq -Q33


# get quality stats
usearch -fastq_stats splitLib/seqs.trimmed.fastq -log splitLib/seqs.stats.log


# remove low quality reads
mkdir qF
usearch -fastq_filter splitLib/seqs.trimmed.fastq -fastq_maxee 0.5 -fastaout qF/seqs.filtered.fasta


# dereplicate sequences. Last step with files separate. 
mkdir deRep
usearch -derep_fulllength qF/seqs.filtered.fasta -output deRep/seqs.filtered.derep.fasta -sizeout


# filter singletons  - This rids sigletons - Decided to do without
# mkdir filterSingles
# usearch -sortbysize deRep/seqs.filtered.derep.fasta -minsize 2 -output filterSingles/seqs.filtered.derep.mc2.fasta


# clusterOTUs
mkdir OTUs
usearch -cluster_otus deRep/seqs.filtered.derep.fasta -otus OTUs/seqs.filtered.derep.repset.fasta


# reference chimera check
mkdir chiCheck
usearch -uchime_ref OTUs/seqs.filtered.derep.repset.fasta -db scripts/gold.fa -strand plus -nonchimeras chiCheck/seqs.filtered.derep.repset.nochimeras.fasta


# label OTUs using puthon script from UPARSE
mkdir labelOTUs
python scripts/fasta_number.py chiCheck/seqs.filtered.derep.repset.nochimeras.fasta OTU_ > labelOTUs/seqs.filtered.derep.repset.nochimeras.otus.fasta


# match original quality filtered reads back to otus - this is with bash derep workaround. 
mkdir matchOTUs
usearch -usearch_global qF/seqs.filtered.fasta -db labelOTUs/seqs.filtered.derep.repset.nochimeras.otus.fasta -strand plus -id 0.97 -uc matchOTUs/otu.map.uc


# make otu table
mkdir otuTable 
python scripts/uc2otutab_mod.py matchOTUs/otu.map.uc > otu-table.txt


# convert to biom  - this step not in aciss, so bring into own computer. 
biom convert --table-type="OTU table" -i otu-table.txt -o otu-table.biom


# **use QIIME 1.7, not 1.8** Dependency problem
# assign taxonomy
assign_taxonomy.py -t gg_13_5_otus/taxonomy/97_otu_taxonomy.txt -r gg_13_5_otus/rep_set/97_otus.fasta -i labelOTUs/seqs.filtered.derep.repset.nochimeras.otus.fasta -o assigned_taxonomy



# add taxonomy to BIOM table
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp assigned_taxonomy/seqs.filtered.derep.repset.nochimeras.otus_tax_assignments.txt -i otu-table.biom -o otu_table.biom


# check sequencing depth.
# print_biom_table_summary.py -i otu_table.biom   ## for qiime <1.8
biom summarize-table -i otu_table.biom -o otu_table_summary.txt  # for qiime >=1.8






