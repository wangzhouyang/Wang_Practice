#!usr/bin/perl -w
use strict;
die "perl $0 [shannon_path_list][out_prefix]" unless @ARGV == 2;
my ($in_f,$out_f) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;
open IN, $in_f or die $!;
open OUT1, ">$out_f.reads_pairs" or die $!;
open OUT2, ">$out_f.base_abundance" or die $!;
open OUT3, ">$out_f.reads_abundance" or die $!;
my %hash;
while(<IN>){
	chomp;
	my $line=$_;
	$line=~/results_10M\/(\S+)\.shannon/;
	my $pre=$1;  
	open INN,"$line"or die $!;
	while(<INN>){
	chomp;
	my @b=split /\t/,$_;
	$hash{$pre}{$b[0]}=$b[1];	
	}
	close INN;
}
close IN;

############################################################################

foreach my $key1 ( sort keys %hash #sort with ASCII order

#		sort with dictionary order
#		sort{(my $da = lc $a) =~ s/[\W_]+//g;
#	 		 (my $db = lc $b) =~ s/[\W_]+//g;
#			 $da cmp $db;
#		}keys %hash

){

	foreach my $key2 (keys %{$hash{$key1}})
	{	
		
        if($key2 eq "reads_pairs"){
            print OUT1 $key1."\t".$hash{$key1}->{$key2}."\n";
		}
		if($key2 eq "base_abundance"){
	        print OUT2 $key1."\t".$hash{$key1}->{$key2}."\n";
	    }
		if($key2 eq "reads_abundance"){
			print OUT3 $key1."\t".$hash{$key1}->{$key2}."\n";
		}
	}
}

close OUT1;
close OUT2;
close OUT3;
