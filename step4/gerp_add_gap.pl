use strict;
use warnings;

open(GFF,$ARGV[0]);
open(LIST,$ARGV[1]);
open(OUT,">$ARGV[2]");
my %hash;
while(<GFF>){
    chomp;
    my @arr = split(/\t/);
    $hash{$arr[3]} = $arr[4];
}
my $i = 1;
while(<LIST>){
    chomp;
    my @arr = split(/\t/);
    if($hash{$i}){
        $i = $hash{$i}+1;
    }
    print OUT $i."\t".$arr[1]."\t".$arr[2]."\n";
    $i++;
}

##perl gerp_add_gap.pl B_ZM4_chr1.gaps.gff3 B_ZM4_9species_use_chr1_gerp_record_add_line_number.txt B_ZM4_chr1_use_gerp.txt