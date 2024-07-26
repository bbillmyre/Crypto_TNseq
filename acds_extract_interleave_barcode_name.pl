use String::Approx 'amatch';
use String::Approx 'aindex';
$first =@ARGV[0];
$second = @ARGV[1];
$out = @ARGV[2];

$acds= "ATAAAACTAACAAAATCGGTTATACGATAACGGTCGGTACGGGATTTTCCCATCCTACTTTCATCCCTG";
$linker= "GTAATAACGACTCACTATAGGGTT...AA....AA...CTCCGCTTAAGGGAC";

`/n/apps/CentOS7/install/pyenv-1.0.0/pyenv/shims/interleave-reads.py $first $second -o temp_inter.fastq`;

$out_1 = $out ."_1.fastq";
$out_2 = $out ."_2.fastq";
$out_3 = $out .".counts.tsv";

open(OUT1,">$out_1");
open(OUT2,">$out_2");

open(INTER,"temp_inter.fastq");

$i=0;
while(<INTER>){
chomp();
if($i%8 == 0){
$name = $_;
}

if($i%8 == 1){
$seq = $_;
}

if($i%8 == 2){
$line_3= $_;
}

if($i%8 == 3){
$qual = $_;
}

if($i%8 == 4){
$sec_name = $_;
}

if($i%8 == 5){
$sec_seq = $_;
}

if($i%8 == 6){
$sec_line_3= $_;
}

if($i%8 == 7){
$sec_qual = $_;

$match= index($seq,$acds);

if($match >-1 && length($seq)>0){


if($sec_seq=~ m/$linker/ && length($sec_seq)>0){

$seq = substr($seq,$match+69);
$qual = substr($qual,$match+69);

$sec_match="@+";
$sec_fir_match="@-";


$barcode=substr($sec_seq,$sec_fir_match,53);
$sec_seq=substr($sec_seq, $sec_match + 2);
$sec_qual=substr($sec_qual, $sec_match + 2);

print OUT1 "\@$barcode\n$seq\n$line_3\n$qual\n";
print OUT2 "\@$barcode\n$sec_seq\n$sec_line_3\n$sec_qual\n";
}

}

}

$i++;

}

`rm temp_inter.fastq`;

$realigned = $out . "picsort.bam";
$counts_out = $out. "_counts.bed";
$bed_sorted = $out. "_sorted.bed";


`perl ../../mapper/assembler_H99.pl $out_1 $out_2 $out`;

`perl count_start_pos.pl $realigned $out > $counts_out`;

`bedtools sort -i $counts_out > $bed_sorted`;

$bed_norm = $out."_normalized.bed";

open(BED,$bed_sorted);
$total_inserts=0;

while(<BED>){


@cur_line = split("\t",$_);


$total_inserts=$total_inserts+@cur_line[4];


}

print "TOTAL INSERTS = $total_inserts\n";

close(BED);

open(BED,$bed_sorted);
open(BED_NORM,">>$bed_norm");

while(<BED>){

@cur_line = split("\t",$_);

$norm_inserts= @cur_line[4]/$total_inserts;

print BED_NORM "@cur_line[0]\t@cur_line[1]\t@cur_line[2]\t$norm_inserts\t@cur_line[4]\t@cur_line[5]" ;

}

`perl add_col.pl $bed_norm Chr_Library $out $out TRUE > $out."col"`;

`perl assign_gene.pl $out."col" crypto_h99_genes_whole_cds.gff $out."labeled.bed"`;







