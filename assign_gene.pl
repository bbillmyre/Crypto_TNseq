

$in = @ARGV[0];
$coords = @ARGV[1];
$out_file = @ARGV[2];
open(OUT, ">>",$out_file);
open(IN, $in);
$i=0;
$total_lines=`wc -l $in`;
chomp($total_lines);
while(<IN>){

	
chomp();
@inserts = split("\t", $_);


#print OUT "TEST\n";

#print "run number $i of $total_lines started\n";

$insert_coord= @inserts[1];
$insert_chr = @inserts[0];

open(COORDS,$coords);

$found=0;

while(<COORDS>){
chomp();
@gene= split("\t",$_);
$gene_chr= @gene[0];
$gene_start=@gene[3];
$gene_end=@gene[4];


if($insert_chr eq $gene_chr && (abs($insert_coord) > $gene_start) && (abs($insert_coord) < $gene_end)){
print "FOUND\t";
print "@gene\t@inserts\n";
$found=1;

$gene_name = @gene[8];

push(@inserts,$gene_name);
print OUT join("\t",@inserts),"\n";
last;
}

}
if(!$found){
print "INTERGENIC\n";
push(@inserts,"Intergenic");
print OUT join("\t",@inserts),"\n";

}
close(COORDS);



}
close(OUT);
