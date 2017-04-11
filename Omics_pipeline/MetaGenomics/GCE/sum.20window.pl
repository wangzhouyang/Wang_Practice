#!/usr/bin/perl

my $usage = "usage: perl $0 targetDirectory sample fileExtension cutoff outputPrefix\n";

die $usage if @ARGV <5;
my ($dir,$sample,$extension,$cutoff,$outPrefix)=@ARGV;

$dir =~ s/\/$//;

my @list = `ls $dir/*.$sample.$extension`;

open STAT,"> $outPrefix.stat" or die "can not access $outPrefix.stat. $!\n";
open OUT,"> $outPrefix.total" or die "can not access $outPrefix.total. $!\n";
my (@reads20,@ref20);
                            
while(@list > 0 ){         
	my $file = shift @list;
	open IN,"< $file" or die "$file can not be access. $!\n";
	my @lines = <IN>;      
#	my $header = shift @lines;   
	my $footer = pop @lines;my @gcsize = split(/\t/,$footer);
	my ($sumReads,$sumRef);
	for(1..20){
		my @cols = split (/\t/,$lines[$_]);
		$sumReads += $cols[2];
		$sumRef += $cols[3];
	}
	next if $sumReads < $cutoff;
	for(1..20){
		my @cols = split (/\t/,$lines[$_]);
		$reads20[$_] += $cols[2];
		$ref20[$_] += $cols[3];
	}
	print STAT "$sumReads\t$sumRef\t$gcsize[2]\t$file";
}
close STAT;

for(1..20){
	my $rate = ($ref20[$_]==0)?0:$reads20[$_]/$ref20[$_]; 
	print OUT "$_\t$reads20[$_]\t$ref20[$_]\t$rate\n";
}

