
use strict;
use warnings;

my $id;
my $new_id;
my $num = 0;

open(LIST,$ARGV[0]);

while(<LIST>){
    chomp;
    if(!/GWAP/){next;}
    my @arr = split(/\t/);
    if($arr[2] eq "gene"){
        $arr[8]=~s/ID=//;
        $arr[8]=~s/;.*//;
        $id = $arr[8];
	$num++;
        $num=sprintf "%06d",$num;
        $new_id = $ARGV[1]."".$num;
    }
    s/$id/$new_id/g;
    print $_."\n";
}
