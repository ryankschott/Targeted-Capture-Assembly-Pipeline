#!/bin/bash

# Targeted capture of complete coding regions across divergent species
# Ryan K Schott, Bhawandeep Panesar, Daren C Card, Matthew Preston, Todd A Castoe, Belinda SW Chang

clean_up() {
	#Deletes broken files if they were created and ends program
	printf "Creating $1" >&2
	rm -f $1
	shift
	for file in $@
	do
		printf ", $file" >&2
		rm -f $file
	done
	printf " failed\n" >&2
	exit 1
}

USAGE="Usage: bash STAMPY_SeqCap_Assembly.sh [-s scriptdir] ref_genome.fas threads [fwd_read.fastq rvrse_read.fastq]..."
scriptdir=$(pwd)
seqs=1

while getopts ':hs:' flag
do
	case $flag in
		h) # Help
			echo "$USAGE"
			exit 1
			;;
		s) # Change script directory
			scriptdir=$OPTARG
			;;
		\?)# Unexpected flag
			echo "Unexpected option -$OPTARG" >&2
			exit 1
			;;
		:) # Missing argument
			echo "Option -$OPTARG requires an argument" >&2
			exit 1
			;;
		*) # Unexpected error
		    echo "Unexpected error in getopts" >&2
			exit 1
			;;
	esac
done
shift $( expr $OPTIND - 1 )

if [ $# -le 1 ]
then
	echo $USAGE
	exit 1
fi

date +"Time started: %F %T"
echo "Using Reference Genome: $1"
ref=$1
shift

echo "Using $1 threads"
threads=$1
shift

if [ -f $ref.stidx ]
then
	echo "Reference exists, skipping"
else
	echo "Building Reference File"
	$scriptdir/stampy.py -G $ref $ref || clean_up $ref.stidx
fi

if [ -f $ref.sthash ]	
then
	echo "Hash table exists, skipping"
else	
	echo "Creating hash table"
	$scriptdir/stampy.py -g $ref -H $ref || clean_up $ref.sthash
fi

if [ -f $ref.fai ]
then
   echo "Samtools index exists, skipping"
else
   echo "Making Samtools Index"
   samtools faidx $ref || clean_up $ref.fai
fi

if [ $# -gt 0 ]
then
	tmp=$( dirname $ref | cut -d / -f 1 )
	dirs=("." "..")
	dirs+=( $(ls -F | grep ".*/" | sed 's|/||') )
	for d in ${dirs[@]}
	do
		if [ $tmp = $d ]
		then
			ref="../../$ref"
			break
		fi
	done

	tmp=$(basename $ref)
	refname=${tmp%.*}
	echo "Making ${refname}_STAMPY directory"
	if [ ! -d ${refname}_STAMPY ]
	then
		mkdir ${refname}_STAMPY
	fi
	cd ${refname}_STAMPY} || { echo "Changing Directory failed"; exit 1; }
fi

while test ${#} -gt 0
do
	date +"Pair start: %T"
	echo "Making Directories"

	tmp=$(basename $1)
	seqname=${tmp%_R*}
	echo "Using Sequence Name: $seqname"
	if [ ! -d "$seqname" ]
	then
		mkdir $seqname
	fi

	forward=$1
	echo "Forward read: $forward"
	shift
	reverse=$1
	echo "Reverse read: $reverse"
	shift
	echo "Changing to directory $seqname"
	cd $seqname || { echo "Changing Directory failed"; exit 1; }

	if [ -f $seqname.sam ]
	then
	   echo "STAMPY output exists, skipping alignment"
	else
	   echo "Running STAMPY"
		$scriptdir/stampy.py -g $ref -h $ref -M $forward $reverse > $seqname.sam || clean_up $seqname.sam
	fi

	then
	   echo "BWA output exists, skipping alignment"
	else
	   echo "Running BWA: bwa mem -B 2 -M -t $threads $ref $forward $reverse"
	   bwa mem -B 2 -M -t $threads $ref $forward $reverse > $seqname.sam || clean_up $seqname.sam
	fi

	if [ -f fixmate_$seqname.bam ]
	then
	   echo "BAM exists, skipping alignment"
	else
	   echo "Fix and Convert to BAM"
	   samtools fixmate -O bam -@ $threads $seqname.sam fixmate_$seqname.bam || clean_up fixmate_$seqname.bam
	fi

	if [ -f sorted_$seqname.bam ]
	then
	   echo "Sorted BAM exists, skipping alignment"
	else
	   echo "Sorting BAM file"
	   samtools sort -@ $threads -o sorted_$seqname.bam fixmate_$seqname.bam || clean_up sorted_$seqname.bam
	fi

	if [ -f sorted_$seqname.bam.bai ]
	then
	    echo "Sorted BAM Index exists, skipping"
	else
	    echo "Indexing Sorted BAM file"
	    samtools index -@ $threads sorted_$seqname.bam || clean_up sorted_$seqname.bam.bai
	fi

	if [ -f stats_$seqname.txt ]
	then
	   echo "Stats outputs exist, skipping"
	else
	   echo "Generating Stats Files"
	   samtools flagstat -@ $threads sorted_$seqname.bam > flagstat_$seqname.txt || clean_up flagstat_$seqname.txt
	   samtools stats -@ $threads sorted_$seqname.bam > stats_$seqname.txt || clean_up stats_$seqname.txt
	   samtools idxstats sorted_$seqname.bam > idxstats_$seqname.txt || clean_up idxstats_$seqname.txt
	   samtools depth -d 1000000 -a sorted_$seqname.bam > depth_$seqname.txt || clean_up depth_$seqname.txt
	fi
	if [ -f depth_of_coverage_$seqname.csv ]
	then
		echo "Depth of Coverage exists, skipping"
	else
		echo "Calculating depth of coverage"
		python $scriptdir/depth_of_coverage_impl.py -b sorted_$seqname.bam -o depth_of_coverage_$seqname.csv -t || clean_up depth_of_coverage_$seqname.csv
	fi

	if [ -f mapped_$seqname.bam ]
	then
	    echo "Mapped output exists, skipping"
	else
	    echo "Filtering out only the mapped seqs"
	    samtools view -b -F 4 -@ $threads sorted_$seqname.bam > mapped_$seqname.bam || clean_up mapped_$seqname.bam
	fi

	if [ -f mapped_$seqname.bam.bai ]
	then
	    echo "Mapped BAM Index exists, skipping"
	else
	    echo "Indexing Mapped BAM file"
	    samtools index -@ $threads sorted_$seqname.bam || clean_up sorted_$seqname.bam.bai
	fi

	if [ -f consensus_$seqname.fq ]
	then
	    echo "consensus fq exists, skipping"
	else
	    echo "Generating consensus with mpileup"
		samtools mpileup -Q 20 -q 5 -d 5000 -uf $ref mapped_$seqname.bam | bcftools call -c | perl $scriptdir/vcfutils.pl vcf2fq -d 5 -Q 20 -l 1 > consensus_$seqname.fq || clean_up consensus_$seqname.fq
	fi

	if [ -f consensus_$seqname.fasta ]
	then
	    echo "consensus fasta exists, skipping"
	else
	    echo "Converting FASTQ to FASTA"
	    python $scriptdir/fastqtofasta.py consensus_$seqname.fq consensus_$seqname.fasta || clean_up consensus_$seqname.fasta
	fi

	if [ -f consensus_uppercase_$seqname.fasta ]
	then
		echo "Uppercase only fasta exists, skipping"
	else
	    echo "Removing Lowercase"
	    perl -e 'while(<>) { if ($_ =~ /^>.*/) { print $_; } else { $_ =~ tr/acgtryswkmbdh/N/; print $_;}}' < consensus_$seqname.fasta > consensus_uppercase_$seqname.fasta || clean_up consensus_uppercase_$seqname.fasta
	fi

	if [ -f consensus_uppercase_blasted_$seqname.xml ]
	then
		echo "Blasted Uppercase consensus exists, skipping"
	else
		echo "Blasting Uppercase consensus file"
		blastn -db /media/Storage2/BlastDB/nt -query consensus_uppercase_$seqname.fasta -outfmt 5 -task dc-megablast -max_target_seqs 1 -num_threads $threads -out consensus_uppercase_blasted_$seqname.xml || clean_up consensus_uppercase_blasted_$seqname.xml
	fi

	if [ -f consensus_uppercase_annotated_$seqname.fasta ]
	then
		echo "Annotated Uppercase consensus exists, skipping"
	else
		echo "Annotating Uppercase consensus file"
		python $scriptdir/seqCaptAnnotateV2.py -b consensus_uppercase_blasted_$seqname.xml -f consensus_uppercase_$seqname.fasta -o consensus_uppercase_annotated_$seqname.fasta || clean_up consensus_uppercase_annotated_$seqname.fasta
	fi

	if [ -f coverage_uppercase_$seqname.csv ]
	then
		echo "Coverage Uppercase file exists, skipping"
	else
		echo "Calculating Uppercase Coverage"
		python $scriptdir/CoveragePercentV2.py -r $ref -i consensus_uppercase_$seqname.fasta -o coverage_uppercase_$seqname.csv || clean_up coverage_uppercase_$seqname.csv
	fi

	if [ -f consensus_ACGT_$seqname.fasta ]
	then
		echo "ACGT only fasta exists, skipping"
	else
	    echo "Removing Lowercase"
	    perl -e 'while(<>) { if ($_ =~ /^>.*/) { print $_; } else { $_ =~ tr/acgtryswkmbdhRYSWKMBDH/N/; print $_;}}' < consensus_$seqname.fasta > consensus_ACGT_$seqname.fasta || clean_up coverage_ACGT_$seqname.fasta
	fi

	if [ -f consensus_ACGT_blasted_$seqname.xml ]
	then
		echo "Blasted ACGT consensus exists, skipping"
	else
		echo "Blasting ACGT consensus file"
		blastn -db /media/Storage2/BlastDB/nt -query consensus_ACGT_$seqname.fasta -outfmt 5 -task dc-megablast -max_target_seqs 1 -num_threads $threads -out consensus_ACGT_blasted_$seqname.xml || clean_up consensus_ACGT_blasted_$seqname.xml
	fi

	if [ -f consensus_ACGT_annotated_$seqname.fasta ]
	then
		echo "Annotated ACGT consensus exists, skipping"
	else
		echo "Annotating ACGT consensus file"
		python $scriptdir/seqCaptAnnotateV2.py -b consensus_ACGT_blasted_$seqname.xml -f consensus_ACGT_$seqname.fasta -o consensus_ACGT_annotated_$seqname.fasta || clean_up consensus_ACGT_annotated_$seqname.fasta
	fi

	if [ -f coverage_ACGT_$seqname.csv ]
	then
		echo "Coverage ACGT file exists, skipping"
	else
		echo "Calculating ACGT Coverage"
		python $scriptdir/CoveragePercentV2.py -r $ref -i consensus_ACGT_$seqname.fasta -o coverage_ACGT_$seqname.csv || clean_up coverage_ACGT_$seqname.csv
	fi

	date +"Pair end: %T"
	echo "_________________________________________________"
	echo "FINISHED. Moving on to next pair"
	echo "_________________________________________________"
	cd ..  || { echo "Changing Directory failed"; exit 1; }
	seqs=$((seqs+1))
done

if [ $seqs -gt 0 ]
then
	echo "Creating Summary Files"
        Rscript $scriptdir/completeness_depth_sumV2.R coverage_uppercase_.*csv Summary_CompletenessUpper.csv coverage_ACGT_.*csv Summary_CompletenessACGT.csv depth_of_coverage_.*csv Summary_Depth.csv || clean_up Summary_CompletenessUpper.csv Summary_CompletenessACGT.csv Summary_Depth.csv

	echo "Summarizing Results"
	bash $scriptdir/Summary_Results.sh || clean_up Summary_Results.csv
fi

echo "All Done here!"
date +"Time ended: %F %T"
