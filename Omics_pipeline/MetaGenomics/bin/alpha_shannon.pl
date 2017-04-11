use warnings;
use strict;

die "perl $0 [*.in] [*.output] [row(optional)]\n" unless @ARGV >= 2;

my ($in_f, $out_f, $row) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;

my (@gene, @sum, @shannon,$title,$shannon) = ();
if($in_f =~/(\/|)(\S+)\.gz$/){
	open IN,"gzip -dc $in_f|" or die $!;
	$title = $2;
}elsif($in_f =~/(\/|)(\S+)$/){
	open IN, $in_f or die $!;
	$title = $2;
}
my @head;
#chomp(my $h = <IN>);
#my @head = split /\s+/, $h;
#shift @head;
while (<IN>) { 
	chomp;
	my @s = split /\s+/;
	if($.==1){
		my $test = (defined $row)?$s[$row]:$s[1];
		if($test =~/\D/ && $test !~ /e-/){
			shift @s;
			@head = @s;
			next;
		}else{
			@head = (1..$#s);
		}
	}
	if (not defined $row){
		shift @s;
		for (0..$#s) {
			next if $s[$_]== 0;
			#$gene[$_]++;
			$sum[$_] += $s[$_];
			$shannon[$_] -= $s[$_] * log($s[$_]);
		}
	}else{
		next if $s[$row]== 0;
		$sum[$row] += $s[$row];
		$shannon[$row] -= $s[$row] * log($s[$row]);
	}
}
close IN;

open OT, ">$out_f" or die "$!\n$out_f";
if(defined $row){
	$shannon[$row] = $shannon[$row] / $sum[$row] +log($sum[$row]);
	print OT "$title\t$shannon[$row]\n";
}else{
	for(0..$#head){
	    $shannon[$_] = $shannon[$_] / $sum[$_] + log($sum[$_]);
		print OT "$head[$_]\t$shannon[$_]\n";
	}
}
close OT;

print STDERR "Program End...\n";
