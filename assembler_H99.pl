
##entry arguments
#################
$pair1 = $ARGV[0]; #first paired end fastq
$pair2 = $ARGV[1]; #second paired end fastq
$nameout = $ARGV[2]; #base name for output/temp files
##################

#paths to ref and jars NEED TO EDIT BASED ON LOCATION/SYSTEM
#################
#$ref = "~/genomes/pombe/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa";
$ref = "/home/rb2091/genomes/cryptococcus/neoformans/FungiDB-44_CneoformansH99_Genome.fasta";
#server picard
#$picard = "/n/local/bin/picard.jar";
##local picard
$picard="/n/apps/CentOS7/bin/picard.jar";
#$picardsort = "~/applications/picard-tools-1.115/SortSam.jar";
#$picardadd = "~/applications/picard-tools-1.115/AddOrReplaceReadGroups.jar";
#$picardmark = "~/applications/picard-tools-1.115/MarkDuplicates.jar";
$gatkjar = "/home/rb2091/apps/GATK/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar";
#################

#naming block for temp files 
#######################
$name1 = $nameout . "_1.sai";
$name2 = $nameout . "_2.sai";
$finalname = $nameout . ".sam";
$finalnamebam = $nameout . ".bam";
$finalnamesort = $nameout . "sort";
$finalnamesortin = $nameout. "sort.bam";
$cleaned = $nameout . "clean.bam";
$named = $nameout. "named.bam";
$cleanreorder = $nameout . "cnr.bam";
$snpcalls = $nameout . "calls.vcf";
$sorted = $nameout . "picsort.bam";
$intervals = $nameout . ".intervals";
$dedup = $nameout . "dups.bam";
$metrics = $nameout . "dupmets"; 
$realigned = $nameout . "realign.bam";
#######################





`bwa mem $ref $pair1 $pair2 -t 3 > $finalname`;
`samtools view -S -b $finalname >$finalnamebam`;

`rm $finalname`;




`java -Xmx6024m -jar $picard CleanSam I=$finalnamebam O=$cleaned`;

`rm $finalnamebam`;

`java -Xmx6024m -jar $picard SortSam I=$cleaned O=$sorted SO=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=./`;

`rm $cleaned`;

#`java -Xmx6024m -jar $picard AddOrReplaceReadGroups INPUT=$sorted OUTPUT=$named RGID=$nameout RGLB=$nameout RGPL=illumina RGPU=1 RGSM=$nameout`;

#`rm $sorted`;

#`java -Xmx6024m -jar $picard MarkDuplicates INPUT=$named OUTPUT=$dedup METRICS_FILE=$metrics`; 
#`rm $named`;

#`samtools index $dedup`;
#`java -Xmx6024m -jar $gatkjar -T RealignerTargetCreator -R $ref -I $dedup -o $intervals`;

#`java -Xmx6024m -jar $gatkjar -T IndelRealigner -R $ref -I $dedup -targetIntervals $intervals -o $realigned`;

#`rm $dedup`;
#`rm $intervals`;
#`rm $metrics`;
#`rm $finalnamesortin`;



`samtools index $sorted`;


 

#default SNP calling block, indels and snps, use outside of script 
#`java -jar $gatkjar -T UnifiedGenotyper -R $ref -I $realigned --sample_ploidy 1 -o $snpcalls -glm BOTH`;
 
 




