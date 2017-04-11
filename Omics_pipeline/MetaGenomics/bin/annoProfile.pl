#! /usr/bin/perl

my $usage ="Here is a script to calculate higher-level's profile via annotation and matched profile.
USAGE:	perl $0 annotateList base-profile output (row)
				annotateList	A file contains a row the same with base-profile and a 2nd row the annotation
				base-profile	The basic profile. Generally refers the gene-level profile
				output		Output file path/name
				(row)		(optional)Specific rows(start from 0) to be calculate, delimited by comma. Default is all.
Author: fangchao\@genomics.cn\n";

use strict;
use warnings;
sub readMethod{$_=shift;($_ =~ /\.gz$/)?"gzip -dc $_|":"< $_"};
sub writeMethod{$_=shift;($_ =~ /\.gz$/)?"|gzip > $_":"> $_"};

die $usage if @ARGV<2;
our ($annotation, $prof, $annoProf, $row) = @ARGV;
our (%ANNO, @cols, @rows);

open AN, readMethod($annotation) or die "can't read $annotation:$!\n";
open IN, readMethod($prof) or die "can't read $prof:$!\n";
open OUT, writeMethod($annoProf) or die "can't write $annoProf. $!\n";

while(<AN>){
	chomp;
	my @a = split /\t/;
	$ANNO{'anno'}{$a[0]} = $a[1];
}
close AN;

#open IN, readMethod($prof) or die "can't read $prof:$!\n";
chomp($_=<IN>);
@cols = split /\t/;
if(defined $row){
	my @a=split(",",$row);
	foreach my $r(@a){push @rows,$r}
}else{
	for(my $i=1;$i<=$#cols;$i++){push @rows,$i}
}

while(<IN>){
	chomp;
	my @a = split /\t/;
	my $anno = $ANNO{'anno'}{$a[0]};
	foreach(@rows){  
		if($anno){
			$ANNO{'val'}{$anno}{$_} += $a[$_];
		}
	}
}
close IN;


#open OUT, writeMethod($annoProf) or die "can't write $annoProf. $!\n";
foreach(@rows){
	print OUT "\t$cols[$_]";
}
print OUT "\n";

foreach my $id (sort keys %{$ANNO{'val'}}){
	print OUT $id;
	foreach my $r(@rows){
		print OUT "\t$ANNO{'val'}{$id}{$r}";
	}
	print OUT "\n";
}
close OUT;





