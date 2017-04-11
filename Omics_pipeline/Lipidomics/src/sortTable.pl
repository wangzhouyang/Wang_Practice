#!/usr/bin/perl -w
use strict;

my $usage = "Usage: perl $0 <input_file> <prefix>\n";
die $usage if @ARGV ne 2;
my ($input,$prf) = @ARGV;

my %HS;
open IN,"<$input" || die $!;
my $head = <IN>;
chomp($head);
my @col=split /\t/,$head;
for (my $i=1;$i<=$#col;$i++){
	my $value = $col[$i];
	$value =~ s/$prf//;
	$HS{$value} = $i;
}

print $col[0];
foreach my $i (sort {$a<=>$b} keys %HS){
	print "\t$col[$HS{$i}]";
}
print "\n";

while(<IN>){
	chomp;
	my @a = split;
	print $a[0];
	foreach my $i (sort {$a<=>$b} keys %HS){
		print "\t$a[$HS{$i}]";
	}
	print "\n";
}
