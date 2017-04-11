#!/usr/bin/perl
use strict;

my $usage = "Description:
	A script calculating Richness(gene count) and index-occurance for profile.
Usage:
	perl $0 <profile>\n";

die $usage if @ARGV < 1;
my $prof = $ARGV[0];
$prof =~ /(\S+).profile$/;
my $richFile = $1.".count";
my $occFile = $1.".cover";

open IN,"< $prof" or die "can not access $prof. $!\n";
open RF,"> $richFile" or die "can not access $richFile. $!\n"; 
open OF,"> $occFile" or die "can not access $occFile. $!\n";

my $head = <IN>; chomp($head);
my @col = split(/\t/,$head);

my %RICH;
while(<IN>){
	chomp;
	my @a = split;
	die "line $.,col number $#a is shorter than head $#col." if $#a != $#col;
	my $element = 0 ;
	for(my $i=1;$i<=$#a;$i++){
		if($a[$i] >0){ 
			$RICH{$i} ++;
			$element ++;
		}
	}
	my $occ = $element / $#a;
	print OF "$a[0]\t$occ\n";
}
close IN;
close OF;
for (my $i=1;$i<=$#col;$i++){
	print RF "$col[$i]\t$RICH{$i}\n";
}
close OF;
