#!usr/bin/perl
use warnings;
use strict;
my $filename = <STDIN>;
open( my $in, '<', $filename ) or die "Cant open $filename: $!"; 
my $line_= $filename;
$line_=readline($filename) if $line_=~/>/;
close $filename;
if($line_ =~ /L|M|F|W|K|Q|E|S|P|V|I|Y|H|R|D/){
print "Protein\n"
}
else {
print "Nucleotide\n"
}
