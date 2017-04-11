#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 <origin_fastQ> <50bp_fastQ>\n" if @ARGV<2;
my($fq1,$out)=@ARGV;

my $openMethod = ($fq1 =~ /gz$/)?"gzip -dc $fq1 |":"$fq1";
my $outMethod = ($out =~ /gz$/)?"| gzip > $out ":"|gzip > $out.gz";

open OUT,$outMethod || die "can't \"$outMethod\". $!\n";
open FQ,$openMethod || die "can't open $fq1. $!\n";
while(<FQ>){
	chomp;
	my $id = $_;
	chomp(my $s = <FQ> );
	chomp(my $d = <FQ> );
	chomp(my $q = <FQ> );
	#do trim
	$s = substr($s,0,50);
	$q = substr($q,0,50);

	print OUT "$id\n$s\n$d\n$q\n";
}

close FQ;
close OUT;
