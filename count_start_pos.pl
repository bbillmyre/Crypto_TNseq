use Bio::Cigar;
#feed in a bam file as first arg and a name for the lib as sec arg
 
$in = @ARGV[0];
$name = @ARGV[1];

`samtools view $in -f 0x63 > tmp_forw.sam`;
`samtools view $in -f 0x53 > tmp_rev.sam`;


open(IN, "tmp_forw.sam");


%insert_sites;
while(<IN>){

@line=split("\t",$_);
$read_name= @line[0];

#print $read_name,"\n";

$combo=@line[2]."_".@line[3]."+";

###saves an array of read names instead of numbers. Then count the items in array in next script when I deal with barcodes.
if($insert_sites{$combo}){
$insert_sites{$combo}=$insert_sites{$combo}.",$read_name";


}else{

$insert_sites{$combo}=$read_name;
}

}
close(IN);

open(REV,"tmp_rev.sam");
while(<REV>){

@line=split("\t",$_);
$read_name=@line[0];
$cigar =Bio::Cigar->new($line[5]);

#print "read_length is ",length(@line[9]),"\n";
#print "cigar_ref is ",$cigar->reference_length,"\n"; 

$rev_pos = @line[3]+ $cigar->reference_length - 1;
$combo=@line[2]."_".$rev_pos."-";


if($insert_sites{$combo}){
$insert_sites{$combo}=$insert_sites{$combo}.",$read_name";

}else{

$insert_sites{$combo}=$read_name;
}

}









@sites= keys %insert_sites;
for(@sites){

$cur_site=$_;

if($insert_sites{$_}){

%uniq_barcode=();
$uniq_num=0;

if(substr($_,-1) eq "+"){
$cur_key=$_;
chop($cur_key);
@location=split("_",$cur_key);
#print $insert_sites{$_};


@barcodes=split(",",$insert_sites{$_});
for(@barcodes){
$cur_bar=$_;
if($uniq_barcode{$cur_bar}){

}else{

$uniq_barcode{$cur_bar}=1;
$uniq_num=$uniq_num+1;
}


}

print @location[0],"\t",@location[1],"\t",@location[1],"\t",$insert_sites{$cur_site},"\t",$uniq_num,"\t","+","\n";
}

if(substr($_,-1) eq "-"){
$cur_key=$_;
chop($cur_key);
@location=split("_",$cur_key);

@barcodes=split(",",$insert_sites{$_});


for(@barcodes){
$cur_bar=$_;
if($uniq_barcode{$cur_bar}){

}else{

$uniq_barcode{$cur_bar}=1;
$uniq_num=$uniq_num+1;
}


}




print @location[0],"\t",@location[1],"\t",@location[1],"\t",$insert_sites{$cur_site},"\t",$uniq_num,"\t","-","\n";
}



}
}




`rm tmp_forw.sam`;
`rm tmp_rev.sam`;
