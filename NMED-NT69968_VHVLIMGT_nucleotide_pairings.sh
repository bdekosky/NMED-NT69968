#!/bin/bash

#Usage command:
#bash NMED-NT69968_VHVLIMGT_nucleotide_pairings.sh EXPT_PREFIX_ 3_Nt_IMGT_OUTPUTFILE

#VH:VL pairing and clustering workflow from IMGT 3_Nt-sequences output files
#3_Nt IMGT tab-delimited files derived from MiSeq VH:VL reads 1 and 2 should be combined into a single file (3_IMGT_OUTPUTFILE above)
#This script outputs several files with different aspects of VH:VL pairing information
#VHVLpairs.txt:  contains all paired CDR3 nt junctions
#completeVHVLpairs: contains all paired CDR3 nt junctions along with V(D)J usage info and CDR3 length
#over1read files (several): contains only those CDRH3:CDRL3 nt pairs observed more than once in the completeVHVLpairs file
#unique_pairs: contains one entry for each cluster seed after 96% clustering of CDRH3 junctions with usearch
#unique_pairs_over1read_noalleles.txt: contains the unique_pairs_over1read entries with V(D)J allele identifications omitted
#The usearch clustering program must be executable via the command "usearch" on line 37 ( http://www.drive5.com/usearch/ , Edgar, R. C. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26, 2460–2461, 2010).  This script has been tested with usearch v5.2.32

#NOTE:  Illumina MiSeq file format sequence headers are required.  If analyzing data downloaded from the NCBI Sequence Read Archive (SRA), the original MiSeq file format and headers must first be reconstructed prior to quality filtering and IMGT analysis.
#The format adjustment is necessary because SRA concatenates the 250bp R1 and R2 reads into a single 500bp sequence, but the script below is written to use the R1/R2 information contained in the MiSeq header which is omitted by the SRA file format.
#
#
#For example, an entry downloaded from the SRA (aside from a '#' beginning each line):
#
#@SRR1056423.266 M01012:72:000000000-A55P3:1:1101:13197:1866 length=500
#TTAGATGGGCCCTGCGATGGGCCCTTGGTGGAGGCTGATGAGACGGTGACCAATGTCCCTTGGCCCCAGATATCCAAAACGCCAACCACGCCCCGTGGCTTTACGCAAGGGCCTTTTTTCCCAAAATCGTAATAATACGACTCCCACTCTTCACATTCGCTTTCTCGCAAATCTCTCGCACAGAAATAGAGAGCCGTGTCTTCCGGTCGCAGGCTGTTCATTTGCAGATGCAGCGTGTTCTTGGAATTGTAGTGCGGCCGCGAGGGTGGGAACAGAGTGACAGTGGGGTTGGCCTTGGGCTGACCTAGGACGGTGACCCTGGTCCCAGTCCCGAAGGTGAGACTGCCACTCCATGTTGAGCAGAAATAATCGGCCTCGTCCCCAGTCTGGAGTCCGGTGATGGCCAGGGTGGCTGCCGTGCCAGACTTGGCGCCAGAGACTTGGGCAGGGATGCCTGAGGGTCGCTTGGTATTTTCATTGATGTGGAGGTTGGCGGCTGT
#+SRR1056423.266 M01012:72:000000000-A55P3:1:1101:13197:1866 length=500
#ABBB3FFFBF4AGGCE2GECFGGGHHHHGGGHGGGHHHHHFH5AEE1GFFHHHHHHFHHHHFGHHH0AEGEGHFFGHFHHEGGDGEF1F?EFE?EEEECHHGFHGG/@GEDFEGHHHHGGFG1FFGHHG?F0<1<<GGHEFFGHHEGG<GFH0DDDDDHGD.EFCGGGC.-CGHFGFE?@BBFFGFFGGG0;AFF?GB..;B;FD?9->--.9AFF./BBFFFFFFFF/BFFFB9;D??BFFFF/;BBB/1>1>11A11>DDA0EEAFGEFCFFH001BEGGF1BFFCE/>EBG1FFHGHGF0BFG1BG1CEE?FABGBFHHHGH01F1B<F<E/</BGBGDGFFHG1CH0F11?22GHHFBCBFCF1?DGFCCAA0.<<<CHH.<<00/<<GF0E@@C.00:C//:GBBA@FFF.;AC?EFGB..;0C0;@B??=-A-B9//9-;B---;9-9/;BFFEBA;9-----//;9/:9/////////9;---;B9/;>@@99
#
#
#The original fastq entries for the above are (aside from a '#' beginning each line):
#
#@M01012:72:000000000-A55P3:1:1101:13197:1866 1:N:0:8
#TTAGATGGGCCCTGCGATGGGCCCTTGGTGGAGGCTGATGAGACGGTGACCAATGTCCCTTGGCCCCAGATATCCAAAACGCCAACCACGCCCCGTGGCTTTACGCAAGGGCCTTTTTTCCCAAAATCGTAATAATACGACTCCCACTCTTCACATTCGCTTTCTCGCAAATCTCTCGCACAGAAATAGAGAGCCGTGTCTTCCGGTCGCAGGCTGTTCATTTGCAGATGCAGCGTGTTCTTGGAATTGT
#+
#ABBB3FFFBF4AGGCE2GECFGGGHHHHGGGHGGGHHHHHFH5AEE1GFFHHHHHHFHHHHFGHHH0AEGEGHFFGHFHHEGGDGEF1F?EFE?EEEECHHGFHGG/@GEDFEGHHHHGGFG1FFGHHG?F0<1<<GGHEFFGHHEGG<GFH0DDDDDHGD.EFCGGGC.-CGHFGFE?@BBFFGFFGGG0;AFF?GB..;B;FD?9->--.9AFF./BBFFFFFFFF/BFFFB9;D??BFFFF/;BBB/
#@M01012:72:000000000-A55P3:1:1101:13197:1866 2:N:0:8
#AGTGCGGCCGCGAGGGTGGGAACAGAGTGACAGTGGGGTTGGCCTTGGGCTGACCTAGGACGGTGACCCTGGTCCCAGTCCCGAAGGTGAGACTGCCACTCCATGTTGAGCAGAAATAATCGGCCTCGTCCCCAGTCTGGAGTCCGGTGATGGCCAGGGTGGCTGCCGTGCCAGACTTGGCGCCAGAGACTTGGGCAGGGATGCCTGAGGGTCGCTTGGTATTTTCATTGATGTGGAGGTTGGCGGCTGT
#+
#1>1>11A11>DDA0EEAFGEFCFFH001BEGGF1BFFCE/>EBG1FFHGHGF0BFG1BG1CEE?FABGBFHHHGH01F1B<F<E/</BGBGDGFFHG1CH0F11?22GHHFBCBFCF1?DGFCCAA0.<<<CHH.<<00/<<GF0E@@C.00:C//:GBBA@FFF.;AC?EFGB..;0C0;@B??=-A-B9//9-;B---;9-9/;BFFEBA;9-----//;9/:9/////////9;---;B9/;>@@99

#Quality filtering of FASTQ data files was performed using the FASTX_toolkit (http://hannonlab.cshl.edu/fastx_toolkit/), with thresholds set at a minimum Q score of 20 over 50 percent of the read.  Quality-filtered data was subsequently converted to FASTA and uploaded to IMGT to obtain the 3_Nt-sequences IMGT files required to run this script.


EXPTNAME=$1
FILE1=$2

rm "$EXPTNAME"unique_pairs.txt

echo "$(date +%b-%d-%H:%M:%S)   Preparing files..."
cat $FILE1 | sed 's/ 1:/_1:/g' | sed 's/ 2:/_2:/g' > pairings_temp_ntdatafile.txt
#print heavy chains in tab-sep file SeqID JunctionNT V-GeneAllele J-GeneAllele D-GeneAllele  CDRH3-length
grep -v '^S' pairings_temp_ntdatafile.txt | awk 'BEGIN {FS="\t"}; {split($4,Vgene," "); split($5,Jgene," "); split($6,Dgene," "); if($3=="productive") if(length($7)>1) print $2 "\t" $16 "\t" Vgene[2] "\t" Jgene[2] "\t" Dgene[2] "\t" length($15)}' | sort -t $'\t' > "$EXPTNAME"IMGT_ntheavys.csv
#same with lights SeqID JunctionNT V-GeneAllele J-GeneAllele CDRL3-length
grep -v '^S' pairings_temp_ntdatafile.txt | awk 'BEGIN {FS="\t"}; {split($4,Vgene," "); split($5,Jgene," "); if($3=="productive") if(length($8)>1) print $2 "\t" $16 "\t" Vgene[2] "\t" Jgene[2] "\t" length($15)}' | awk 'BEGIN {FS="_"}; {if (substr($2,1,1)==1) print $1 "_" "2" substr($2,2); else if (substr($2,1,1)==2) print $1 "_" "1" substr($2,2)}' | sort -t $'\t' > "$EXPTNAME"IMGT_ntlights_readrev.csv

#join for VH-VL junction pairings file
join -t $'\t' -j 1 -o 1.2 2.2 "$EXPTNAME"IMGT_ntheavys.csv "$EXPTNAME"IMGT_ntlights_readrev.csv | sort | uniq -c | sort -n -r | awk '{if ($2!=$3) print $0}' | sed "s/^ *//" | sed 's/ /\t/g' > "$EXPTNAME"VHVLpairs.txt
#join for complete info pairings file HjunctionNT LjunctionNT HV-gene HJ-gene HD-gene HCDR3-length LV-gene LJ-gene LCDR3-length
join -t $'\t' -j 1 -o 1.2 2.2 1.3 1.4 1.5 1.6 2.3 2.4 2.5 "$EXPTNAME"IMGT_ntheavys.csv "$EXPTNAME"IMGT_ntlights_readrev.csv | sort | uniq -c | sort -n -r | awk '{if ($2!=$3) print $0}' | sed "s/^ *//" | sed 's/ /\t/g' > "$EXPTNAME"completeVHVLpairs.txt

#echo "$(date +%b-%d-%H:%M:%S)   Creating a list of unique paired nt junctions read more than 1 times..."
awk -v exptname="$EXPTNAME" 'BEGIN {i=1}; {if ($1>1) print ">" exptname "_completepairs_CDR-H3_Rank_" i "_\n" $2; i=i+1}' "$EXPTNAME"VHVLpairs.txt > "$EXPTNAME"heavynt_junctions.fasta

#echo "$(date +%b-%d-%H:%M:%S)   96% clustering the CDR-3 junction sequences..."
usearch -cluster "$EXPTNAME"heavynt_junctions.fasta -w 4 --maxrejects 0 -usersort --iddef 2 --nofastalign -id 0.96 -minlen 11 -uc results.uc -seedsout seeds.uc

#generate list of clonotypes & frequencies based on unique sequences
awk '/>/{getline;print}' seeds.uc > "$EXPTNAME"uniqueHCntseqs.txt

while read line
do      
    grep -m 1 "$line" "$EXPTNAME"completeVHVLpairs.txt >> "$EXPTNAME"unique_pairs.txt
done < "$EXPTNAME"uniqueHCntseqs.txt

#Filter unique_pairs for seqs with >=2 reads
awk '{if($1>1)print}' "$EXPTNAME"unique_pairs.txt | sort -n -r > "$EXPTNAME"unique_pairs_over1read.txt

#Generate noallele unique_pairs files for Excel pivot table heat maps
awk '{gsub("\*[0-9]*\t", "\t"); print}' "$EXPTNAME"unique_pairs_over1read.txt > "$EXPTNAME"unique_pairs_over1read_noalleles.txt

#Output FASTA files with CDR-H3 for unique pairs and top heavy/light filtered pairs
awk -v exptname="$EXPTNAME" 'BEGIN {i=1}; {print ">" exptname "uniquepairs_CDR-H3_Rank_" i "_\n" $2; i=i+1}' "$EXPTNAME"unique_pairs.txt > "$EXPTNAME"CDR-H3_unique_seqs.fasta
awk -v exptname="$EXPTNAME" 'BEGIN {i=1}; {print ">" exptname "uniquepairs_CDR-H3_Rank_" i "_\n" $2; i=i+1}' "$EXPTNAME"unique_pairs_over1read.txt > "$EXPTNAME"CDR-H3_unique_seqs_over1read.fasta

echo "$(date +%b-%d-%H:%M:%S)   Job complete."

echo 'head ' "$EXPTNAME"'VHVLpairs.txt'
head "$EXPTNAME"VHVLpairs.txt
echo 'Total number of paired reads:'
awk '{sum+=$1}; END {print sum}' "$EXPTNAME"completeVHVLpairs.txt
echo "$EXPTNAME"'completeVHVLpairs.txt'

rm pairings_temp_ntdatafile.txt
rm "$EXPTNAME"IMGT_ntheavys.csv
rm "$EXPTNAME"IMGT_ntlights_readrev.csv

wc -l "$EXPTNAME"unique_pairs.txt
wc -l "$EXPTNAME"unique_pairs_over1read.txt
