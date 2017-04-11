#!/usr/bin/perl -w
use strict;
use POSIX;

die "Usage: perl $0 <gene abundance> <mgs cluster> <output>\n" if @ARGV < 3;
my($abun,$mgs,$out) = @ARGV;
my %HASH;

open AB,"<$abun" or die "cannot open $abun. $!\n";
open CL,"<$mgs" or die "cannot open $mgs. $!\n";
open OT,">$out" or die "cannot open $out. $!\n";

while(<AB>){
	chomp;
	my @a = split;
	$HASH{$a[0]} = $a[-1];
}
close AB;
while(<CL>){
	chomp;
	my @a = split /\t|\s/;
	my $id=shift @a; 
	my $len=scalar(@a);
	my @arr;
	foreach my $g(@a){ push @arr,$HASH{$g}}
	my $val = &q3(@arr);
	print OT "$id\t$val\n";
}
close CL;
close OT;

#####################
#####################
sub q3 {
	my @array = sort {$a<=>$b} @_;
	my $len = scalar(@array);
	my $n3 = $len - int($len/2 + 1.5)/2;
	my $q3 = ($array[floor($n3)] +  $array[ceil($n3)])/2;
	return($q3);
}


