# Targeted-Capture-Assembly-Pipeline

Targeted capture of complete coding regions across divergent species

Ryan K Schott, Bhawandeep Panesar, Daren C Card, Matthew Preston, Todd A Castoe, Belinda SW Chang

NOTE: This is an updated version of the BWA and STAMPY assembly pipelines from the paper

Instructions for Assembly

Prerequisites: Linux machine with BWA and NGM installed and accessible from the path for their respective pipelines. Samtools also needs to be installed and accessible from the path. STAMPY should be present in the current directory. The BLASTn nucleotide database will need to be downloaded to perform BLAST annotation. Within the pipeline script (eg., BWA.sh) find and replace all '/media/Storage2/BlastDB/nt' with the directory and filename used. All other bundled scripts should be present in the current directory. A reference fasta file for assembly will need to be provided. For text matching to work correctly the fasta headers in the reference file should consist of an unique genus and/or species identifier followed by an underscore and then the gene code (eg., '>genus_gene' or '>species_gene', but not '>genus_species_gene')

1. Trim all reads using for example Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)

2. Select either the BWA or STAMPY1 pipeline

3. Run the pipeline:

./BWA_SeqCap_Assembly_v2.sh <Reference Fasta File> <Number of CPU cores to use> <Forward read(s) sp1> <Reverse read(s) sp1> <Forward read(s) sp2> <Reverse read(s) sp2> ... &> <Logfile> & disown -a

Wildcard characters are supported and can be used to run multiple species sequentially as long as the forward and reverse reads names match up (eg., Species1_R1.fastq, Species1_R2.fastq, Species2_R1.fastq, Species2_R2.fastq...).

e.g., ./BWA_SeqCap_Assembly_v2.sh Anolis_reference.fas 20 /path/to/reads/*.fastq &> Logfile.log & disown -a
