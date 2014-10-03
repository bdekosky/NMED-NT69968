#!/bin/bash

#Usage command:
#bash NMED-NT69968_allelic_inclusion_finder.sh replicate1_unique_pairs_over1read_noalleles_file replicate2_unique_pairs_over1read_noalleles_file replicate1_complete_pairs_file replicate2_complete_pairs_file

#This script uses the file output from NMED-NT69968_VHVLIMGT_nucleotide_pairings.sh (specifically, complete_pairs.txt and unique_pairs_over1read_noalleles.txt) to detect heavy chains that paired with multiple light chains, where the allelically included pairings were detected in both technical replicates.

UNIQUE_FILE1=$1
UNIQUE_FILE2=$2
COMPLETE_FILE1=$3
COMPLETE_FILE2=$4

echo "$(date +%b-%d-%H:%M:%S)   Preparing files..."

#Generate a temporary list of all HC's and LC's of interest from unique_files
awk '{print $2 "\t" $1}' "$UNIQUE_FILE1" | sort | uniq | sort > tempHC1.txt
awk '{print $2 "\t" $1}' "$UNIQUE_FILE2" | sort | uniq | sort > tempHC2.txt
awk '{print $3}' "$UNIQUE_FILE1" | sort | uniq -c | sort -n -r | awk '{if($1>1)print}' > tempLC1.txt
awk '{print $3}' "$UNIQUE_FILE2" | sort | uniq -c | sort -n -r | awk '{if($1>1)print}' > tempLC2.txt

awk '{if($1>1)print}' "$COMPLETE_FILE1" > temp_completefile1.txt
awk '{if($1>1)print}' "$COMPLETE_FILE2" > temp_completefile2.txt

join -j 1 -o 1.1 tempHC1.txt tempHC2.txt > tempHC_1and2.txt

rm 1HC_multipleLC_resultsfile.txt

echo "$(date +%b-%d-%H:%M:%S)   Searching for HC's with multiple LC's in both sets..."

rm counterfile.txt

while read line
do 
	rm tempLClist1.txt
	rm tempLClist2.txt
	grep $line temp_completefile1.txt | awk '{print $3 "\t" $1 "\t" $8 "\t" $9 "\t" $10}' | sort | uniq > tempLClist1.txt
	grep $line temp_completefile2.txt | awk '{print $3 "\t" $1}' | sort | uniq > tempLClist2.txt
	join -j 1 -o 1.2 2.2 1.1 1.3 1.4 1.5 tempLClist1.txt tempLClist2.txt > temp_overlappingLCs.txt
	
	#Iterate if HC paired with more than 3 LC
	HCpaircount1=$(grep $line temp_completefile1.txt | awk '{print $3}' | sort | uniq | wc -l | awk '{print $1}')
	if [ $HCpaircount1 -gt 4 ]; then
	continue
	fi
	HCpaircount2=$(grep $line temp_completefile2.txt | awk '{print $3}' | sort | uniq | wc -l | awk '{print $1}')
	if [ $HCpaircount2 -gt 4 ]; then
	continue
	fi

	#Iterate if only one LC paired with that heavy
	overlappingLCcount=$(wc -l temp_overlappingLCs.txt | awk '{print $1}')
	if [ $overlappingLCcount == 1 ]; then
	continue
	fi

	#Iterate if more than 3 LC seqs paired with that heavy to prevent promiscuously paired HCs
	if [ $overlappingLCcount -gt 3 ]; then
	continue
	fi

	#Iterate if nothing in the file
	if [ $overlappingLCcount == 0 ]; then
	continue
	fi

	#Iterate if only one VL seq paired with that heavy
	numVLseqs=$(awk '{print $3}' temp_overlappingLCs.txt | sort | uniq | wc -l | awk '{print $1}')
	if [ $numVLseqs == 1 ]; then
	continue
	fi

	#Iterate if only one VL gene paired with that heavy
	numVLgenes=$(awk '{print $4}' temp_overlappingLCs.txt | sort | uniq | wc -l | awk '{print $1}')
	if [ $numVLgenes == 1 ]; then
	continue
	fi

	#Iterate if only one VJ gene paired with that heavy
	numVJgenes=$(awk '{print $5}' temp_overlappingLCs.txt | sort | uniq | wc -l | awk '{print $1}')
	if [ $numVJgenes == 1 ]; then
	continue
	fi

	hcseq=$(echo "$line" | awk '{print $1}')
	hccount1=$(echo "$line" | awk '{print $2}')
	hccount2=$(echo "$line" | awk '{print $3}')
	awk -v hcseq="$hcseq" -v count1="$hccount1" -v count2="$hccount2" 'BEGIN {print "HC:  " hcseq "\t" count1 "\t" count2 "\nLCs:"}; {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}; END {print "\n\n"}' temp_overlappingLCs.txt >> 1HC_multipleLC_resultsfile.txt
	
	echo >> counterfile.txt
	COUNTER=$(wc -l counterfile.txt | awk '{print $1}')
	echo "loop $COUNTER"
done < tempHC_1and2.txt

rm counterfile.txt

echo "$(date +%b-%d-%H:%M:%S)   Job complete.  1HCmultipleLC seqs are located in 1HC_multipleLC_resultsfile.txt"
