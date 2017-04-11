#! /bin/usr/perl -w
use strict;

my $usage ="Usage:
	perl <kegg.catalog.with.gene.index.file>  <gene.prof> <output>\n";

die "$usage.$!\n" if @ARGV <3;
my($ko,$index,%HASH,%ABUN);

sub openM{my $f=shift;return(($f=~/\.gz$/)?"gzip -dc $f|":"$f")};
open KO,&openM($ARGV[0]) or die "can not access $ARGV[0]. $!\n";
open AB,&openM($ARGV[1]) or die "can not access $ARGV[1]. $!\n";
open OUT,"> $ARGV[2]" or die "can not write $ARGV[3]. $!\n";
my $heads = <AB>;
print OUT $heads;

while(<KO>){
	chomp;
	my @ko=split;
	chomp($_=<AB>);
	my @ab=split;
	next if @ko < 2 ;
	die "KO:$ko[0] ne AB:$ab[0]" if $ko[0] ne $ab[0];
	for(my $i=1;$i<@ko;$i++){
		for(my $s=1;$s<@ab;$s++){
			$ABUN{$ko[$i]}{$s} += $ab[$s];}
	}
}
close KO;
close AB;

foreach my $k(sort keys %ABUN){
	print OUT $k;
	foreach my $s(sort {$a<=>$b} keys %{$ABUN{$k}}){
		print OUT "\t$ABUN{$k}{$s}"
	}
	print OUT "\n";
}
close OUT;
