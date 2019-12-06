#This script takes paired end read (adaptor trimmed) fastq files, merge two fastq files based on overlap (complementary) of tags, collapse the identical reads and map the reads on the expected amplicon and identify the deleted regions in the read and finally report the micro homolay around the repaired sites.
#Arg1 is the fastq file 1 (adaptor should be removed) "pX330_A_R1_trimmed.fastq"
#Arg2 is the fastq file 2 (adaptor should be removed) "pX330_A_R2_trimmed.fastq"
#Arg3 is sample name "pX330"
#Arg4 is the name of file that has expected amplicon sequence
#Usage CRISPR-DNA-Repair.sh pX330_A_R1_trimmed.fastq pX330_A_R2_trimmed.fastq pX330 subject.onlyseq.fa 

#removing the output file if it exists.
rm $3-align-repeat-motif.txt

#Step1
#joining overlapping reads. The pipeline uses a program called "flash" to do this.

flash $1 $2 -o pX330

#Step2
#Extracting the first and last 10 base from amplicon. After a cut if the cutpoint is repaired then it is expected to have sequences from the begining and end of the amplicon. If the begining or end sequences are missing then that would indicate there was not repair. This step is done to avaoid analyzing those reads that are not useful to study DNA repair. 
first10=`cat $4 | awk '{print substr(\$1,1,10)} '`
last10=`cat $4 | awk '{print substr(\$1,length-9)} '`
SEQ=`cat $4 | awk '{print \$1}'`
SEQLEN=`cat $4 | awk '{print length (\$1)}'`

#Step3
#Counting the perfectly repaired sequence
WTCF=`awk 'NR%4==2' $3.extendedFrags.fastq | grep -w -f $4 | wc -l | awk '{print \$1}'`
echo $WTCF $SEQ WT WT WT WT   | awk '{printf ("%6d\t%152s\t%4s\t%15s\t%4s\t%s\n",$1,$2,$3,$4,$5,$6)}'

#Step4
#Collapsing the identical reads. Also only those reads are being considered that have deletions (no insersions). If the length of read is longer than perfectly repaired (sequence identical to WT) clone then this would indicate while repair additional sequences are added. Also considering only those reads that have 10 bases from from the 5' end 10 bases from 3' end exactly identical to expected amplicon.
echo "awk 'NR%4==2' $1 | awk 'length(\$1)<$SEQLEN' | grep ^$first10 | grep $last10$ | sort | uniq -c | awk 'NF>=2'" > Foo.sh
sh Foo.sh > TeMp.fa

#Step5
#Making subject fasta file that would be used for mapping
cat $4 | awk '{printf (">subject\n%s\n",$1)}' > subject.fa

#Step6
#Reading each collapse read and doing global alignment. To do the alignment pipeline uses needle program from EMBOSS package.
cat TeMp.fa | while read line
	do
	echo "$line" | awk '{printf (">%d\n%s\n",$1,$2)}' > query.fa
	CF=`echo "$line" | awk '{printf ("%d\n", \$1)}'` 
	needle -auto -asequence query.fa -bsequence subject.fa -outfile align1 -gapopen 10.0 -gapextend 0.0 -sid1 UNMAPPED -sid2 DATABASE -awidth3 500

	align_length=`cat align1 | awk '\$1~/UNMAPPED/ && NF==4 {print length(\$3)}'`
	indelstart=`cat align1 | awk '$1~/UNMAPPED/ && NF==4 {print \$3}' | awk 'BEGIN {FS = "-"} {print (length(\$1))+1}'`
	alignedseq=`cat align1 | awk '$1~/UNMAPPED/ && NF==4 {print \$3}'`
	alignpart1=`cat align1 | awk '$1~/UNMAPPED/ && NF==4 {print \$3}' | awk 'BEGIN {FS = "-"} {print length(\$1)}'`
	alignpart2=`cat align1 | awk '$1~/UNMAPPED/ && NF==4 {print \$3}' | awk 'BEGIN {FS = "-"} {print length($NF)}'`
	alignpart2seq=`cat align1 | awk '\$1~/UNMAPPED/ && NF==4 {print \$3}' | awk 'BEGIN {FS = "-"} {print $NF}'`
	afterindelends=`echo "$align_length $alignpart2" | awk '{print (\$1-\$2)+1}'`
	NUMBEROFSPLTINDEL=`echo $alignedseq | sed -e s/-/\ /g | awk '{print NF}'`
	NUMBEROFINDEL=`expr $align_length - $alignpart1 - $alignpart2 `
#Step7
#Considering only those alignment where all the deleted bases are contiguous. 
	if [ $align_length -eq $SEQLEN ] && [ $NUMBEROFSPLTINDEL -eq 2 ];
		then
		for (( j = 15 ; j >=1 ; j-- ))
        	do
		echo "cat $4 | awk '{print substr(\$1,$indelstart,$j)}'" > Foo1.sh
        	START_MOTIF=`sh Foo1.sh`
		echo "cat $4 | awk '{print substr(\$1,$afterindelends,$j)}' " > Foo2.sh
        	END_MOTIF=`sh Foo2.sh`
#Step8
#Matching the motif
        	if [ $START_MOTIF = $END_MOTIF ];
			then
                	echo "$CF $alignedseq $align_length $START_MOTIF $NUMBEROFINDEL $indelstart" | awk '{printf ("%6d\t%152s\t%4d\t%15s\t%4d\t%4d\n",$1,$2,$3,$4,$5,$6)}' | awk '{if ($NF<=4) {printf ("%s\tC-NHEJ\n",$0)} else if ($NF=="WT") {printf ("%s\tWT\n",$0)} else {printf ("%s\tSSA\n",$0)}}' >> $3-align-repeat-motif.txt
                	echo "$CF $alignedseq $align_length $START_MOTIF $NUMBEROFINDEL $indelstart" | awk '{printf ("%6d\t%152s\t%4d\t%15s\t%4d\t%4d\n",$1,$2,$3,$4,$5,$6)}' | awk '{if ($NF<=4) {printf ("%s\tC-NHEJ\n",$0)} else if ($NF=="WT") {printf ("%s\tWT\n",$0)} else {printf ("%s\tSSA\n",$0)}}' 
                	break;
        	elif [ $j == "1" ];
			then
                	echo "$CF $alignedseq $align_length NOMOTIF $NUMBEROFINDEL $indelstart" | awk '{printf ("%6d\t%152s\t%4d\t%15s\t%4d\t%4d\n",$1,$2,$3,$4,$5,$6)}' | awk '{if ($NF<=4) {printf ("%s\tC-NHEJ\n",$0)} else if ($NF=="WT") {printf ("%s\tWT\n",$0)} else {printf ("%s\tSSA\n",$0)}}' >> $3-align-repeat-motif.txt
                	echo "$CF $alignedseq $align_length NOMOTIF $NUMBEROFINDEL $indelstart" | awk '{printf ("%6d\t%152s\t%4d\t%15s\t%4d\t%4d\n",$1,$2,$3,$4,$5,$6)}' | awk '{if ($NF<=4) {printf ("%s\tC-NHEJ\n",$0)} else if ($NF=="WT") {printf ("%s\tWT\n",$0)} else {printf ("%s\tSSA\n",$0)}}' 
                	break;
		fi
	done
	fi
done
#Step9
#Writing the output
awk '$4!="NOMOTIF"' $3-align-repeat-motif.txt | awk '{for (i=1; i<=$1; i=i+1) {print $4}}' | sort | uniq -c | sort -nk1 > $3-motif-freq.txt
