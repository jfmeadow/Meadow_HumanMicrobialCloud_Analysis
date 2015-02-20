#!/bin/sh

#-----------------------------
# R1 index read has first 8bp barcodes
# R2 index read has second 8bp barcodes
# Combine these and paste them onto R1

# R1 index: r1index.fastq
# R2 index: r2index.fastq
# R1 read: r1read.fastq
# R2 read: r2read.fastq

# Put all 4 in same folder called `rawData`
# Navigate to this script in that folder. 
# ./manipulatePickle.sh
#-----------------------------


#-----------------------------
# These files came from DFCI:
# 060314_PB2_AA1064_NoIndex_L001_R1_001.fastq.gz
# 060314_PB2_AA1064_NoIndex_L001_R2_001.fastq.gz
# 060314_PB2_AA1064_NoIndex_L001_R3_001.fastq.gz
# 060314_PB2_AA1064_NoIndex_L001_R4_001.fastq.gz
# 
# and were renamed like this: 
# r1read.fastq
# r1index.fastq
# r2index.fastq
# r2read.fastq
# 
#-----------------------------



#-----------------------------
# For quick run, spit out top of file:
# head -n 4000000 r1index.fastq > ../r1index.fastq
# head -n 4000000 r2index.fastq > ../r1index.fastq
# head -n 4000000 r1read.fastq > ../r1read.fastq
# head -n 4000000 r2read.fastq > ../r2read.fastq

#-----------------------------







#-----------------------------
# combine R1 index and R2 index 
#
# revcomp R2 read  - MAYBE NOT
# fastx_reverse_complement -i r2read.fastq -o r2readRevcomp.fastq -Q33
# echo '\nall done revComping R2\n'
#
# paste R1 and revcomp R2
# paste -d '\0' r1index.fastq r2indexrevcomp.fastq > barcodes.fastq
paste -d '\0' r1index.fastq r2index.fastq > barcodes.fastq
echo '\nall done pasting barcodes\n'
#
# renames R2 reads 
# fastx_renamer -i r2readRevcomp.fastq -o r2readRevcompRenamed.fastq -n COUNT -Q33
# echo "\nall done renaming R2\n"
#
# remove ++
perl -i.bak -p -e 's/^\+\+\n/\+\n/' barcodes.fastq
echo '\nall done cleaning up ++\n'
#
# renames sequences
fastx_renamer -i barcodes.fastq -o barcodesRenamed.fastq -n COUNT -Q33
# echo "\nall done renaming. Next to Doug's script\n"
#
# renames R1 reads 
fastx_renamer -i r1read.fastq -o r1readRENAMED.fastq -n COUNT -Q33
#
# renames R2 reads 
fastx_renamer -i r2read.fastq -o r2readRENAMED.fastq -n COUNT -Q33
#

#-----------------------------






#-----------------------------
# Doug Turnbull's perl script to remove spacers:
#
# Spacer sequences
#   The spacer is actually the first 4 of gene primer for first (0 length spacer)
# These are r1read.fastq
# 319F:
#    sequencingPrimer             spacer  trickPrimer  keepPrimer 
# 5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG ACTC CTAC GGGAGGCAGCAG
# 5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TACTC CTAC GGGAGGCAGCAG
# 5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CTACTC CTAC GGGAGGCAGCAG
# 5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG GGTACTC CTAC GGGAGGCAGCAG
# 5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG AACGACTC CTAC GGGAGGCAGCAG
# 5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TTGTTACTC CTAC GGGAGGCAGCAG

# These are r2read.fastq
# 806R:
#    sequencingPrimer             spacer  trickPrimer  keepPrimer 
# 5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GGAC TAC HVGGGTWTCTAAT
# 5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG AGGAC TAC HVGGGTWTCTAAT
# 5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG TAGGAC TAC HVGGGTWTCTAAT
# 5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG CCCGGAC TAC HVGGGTWTCTAAT
# 5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ATTTGGAC TAC HVGGGTWTCTAAT
# 5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GCACAGGAC TAC HVGGGTWTCTAAT

# phased_read_trimmer_PE.pl <inputfile.fastq> <variable trim seq 1> <variable trim seq 2> <variable trim seq 3> <variable trim seq 4> <variable trim seq 5> <variable trim seq 6> <invariant primer sequence> > <outputfile.fastq>
#

#### Must run this once for each gene 

# R1 = 319F primer
perl ../phase6_read_trimmer_PE.pl r1readRENAMED.fastq ACTC TACTC CTACTC GGTACTC AACGACTC TTGTTACTC CTAC > r1readTRIM.fastq  # make this the gene name - one for each gene. 
echo '\nall done with R2\n'
# R2 = 806R primer
perl ../phase6_read_trimmer_PE.pl r2readRENAMED.fastq GGAC AGGAC TAGGAC CCCGGAC ATTTGGAC GCACAGGAC TAC > r2readTRIM.fastq
echo '\nall done with R1\n'
#

# This script writes even sequnces that fail. 

#-----------------------------


#-----------------------------
# move these files to the QIIME directory
#  Actually the processPickle.sh sources from this dir so not necessary. 
# cp r1readTRIM.fastq r2readTRIM.fastq barcodesRenamed.fastq ../QIIME/
# cp r1readTRIM.fastq barcodesRenamed.fastq ../QIIME/






#-----------------------------
# Join paired ends - also once for each of the different genes. 
# download missing dependencies. 
# http://www.wernerlab.org/software/macqiime/add-fastq-join-to-macqiime-1-8-0
join_paired_ends.py -f r1readTRIM.fastq -r r2readTRIM.fastq -b barcodesRenamed.fastq -o joined/
#-----------------------------

















