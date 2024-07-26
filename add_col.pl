


$in= @ARGV[0];

$header = @ARGV[1];
$entry = @ARGV[2];
$library = @ARGV[3];
$add_header=@ARGV[4];

$library = "\"".$library."\"";
$i=0;
open(IN,$in);

while(<IN>){
chomp();
@line= split("\t",$_);

push(@line,$library);

$line=join("\t",@line);

print $line ,"\n"

}



