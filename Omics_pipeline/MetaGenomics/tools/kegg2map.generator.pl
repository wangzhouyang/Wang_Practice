#!/usr/bin/perl -w
use strict;

my $usage = "Usage:
	perl $0 [map file] [val-index prof] > [key-index prof] 2> [log file]\n";
die "$usage" if @ARGV < 2;
my ($map,$prof) = @ARGV;
my %HASH;
sub openM{$_=shift;return(($_=~/\.gz$/)?"gzip -dc $_":"$_")}

open MAP,&openM($map) or die "can not access $map. $!\n";
open PROF,&openM($prof) or die "can not access $prof. $!\n";

$_ = <PROF>;
@_ = split /\t|\s/;
my $size = $#_;
print $_;
while(<PROF>){
	chomp;
	my @a=split;
	die "line $. find $#a values, which less than head $size." if $#a ne $size;
	for(my $i=1;$i<@a;$i++){
		$HASH{$a[0]}{$i} = $a[$i];
	}
}
close PROF;

while(<MAP>){
	chomp;
	my @a=split;
	my @keys=split(',',$a[2]);
	my (%ABUN) = ();
	my ($size,$hit) = (0,0);
	for(my $k=0;$k<@keys;$k++){
		$size ++;
		next unless $HASH{$keys[$k]};
		$hit ++;
		foreach my $i(sort {$a<=>$b} keys %{$HASH{$keys[$k]}}){
			$ABUN{$i} += $HASH{$keys[$k]}{$i}
		}
	}
	print $a[0];
	print STDERR "$a[0]\t$size\t".$hit/$size."\n";
	if(%ABUN){
		foreach my $i(sort {$a<=>$b} keys %ABUN){
			print "\t$ABUN{$i}";
		}
	}else{print "\t0" x $size}; 

	print "\n";
}
close MAP;
exit;
