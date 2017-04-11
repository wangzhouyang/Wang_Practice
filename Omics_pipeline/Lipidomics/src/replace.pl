#!/usr/bin/perl -w
use strict;

my $usage = " Usage: perl $0 <ID mapping file> <input file > [mode]
		mode	h	<only replace the header>
			r	<only replace the row names>
			b	<replace both column and row>
			a	<(default) scan all the file for replacing>
";

if(@ARGV<2){
	die "insufficient ARGV\n$usage";
}
my ($id_f,$in_f,$mode) = @ARGV;
$mode||="a";

my (%MAP);
open ID,"<$id_f"||die "$!\ncan't open $id_f\n";
while (<ID>){
	chomp;
	my @a=split;
	$MAP{$a[1]}=$a[0];
}
close ID;

open FI,"$in_f"||die "$!\ncan't open $in_f\n";
my $head=<FI>;
my $mapping;
chomp($head);
$head =~ s///;
my @a = split(/\t/,$head);
if(not defined $MAP{$a[0]}){$mapping=$a[0]}else{$mapping=$MAP{$a[0]}}
print $mapping."\t";
for (my $i=1;$i<=$#a;$i++){
	if(not defined $MAP{$a[$i]}){
		$mapping=$a[$i]
	}else{$mapping=$MAP{$a[$i]}}
	if ($mode =~ /[hba]/){
		print $mapping."\t";
	}else{
		print $a[$i]."\t";
	}
}
print "\n";

while (<FI>){
	chomp;
	$_ =~ s/^M//;
	if ($mode eq "h"){
		print $_."\n";
		next;
	}else{
		chomp;
		my @a=split;
		print $MAP{$a[0]}."\t";
		for (my $i=1;$i<=$#a;$i++){
			if ($mode =~/[a]/){
				print $MAP{$a[$i]}."\t";
			}else{
				print $a[$i]."\t";
			}
		}
		print "\n";
	}
}
close FI;

