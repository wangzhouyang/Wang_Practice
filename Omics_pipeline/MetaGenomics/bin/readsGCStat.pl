#!/usr/bin/perl -w
use warnings;
use strict;
use File::Basename; 

die &usage if @ARGV !=2;
sub usage {
	print <<USAGE;
usage:
	perl $0 fq <prefix>
USAGE
	exit;
}

my ($fq,$pfx) = @ARGV;
open FQ,"gzip -dc $fq |" or die "error\n";
open OUT,"|gzip >$pfx.GC.gz" or die "error\n";

my $lineNum = 0;
while(<FQ>){
	chomp;
	my @a =split;
	$_ =~ s/@//;
	chomp(my $seq=<FQ>);
	chomp(my $num=<FQ>);
	chomp(my $quality=<FQ>);
	my $GC_count=0;
	my $N_count=0;
	my $len = length($seq);
	# GC cal
	$GC_count=$seq =~ tr/([GC])/\1/;
	$N_count=$seq =~ tr/N/N/;

	print OUT "$_\t$len\t$GC_count\t$N_count\n";
	# stat
	$lineNum ++;
}
close OUT;

print STDERR "$fq done!\t$lineNum\tlines calculated.\n";

