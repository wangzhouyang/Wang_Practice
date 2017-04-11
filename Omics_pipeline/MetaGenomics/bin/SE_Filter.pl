#!/usr/bin/perl -w
use warnings;
use strict;
use File::Basename; 

die &usage if @ARGV < 7;
sub usage {
	print <<USAGE;
usage:
se pattern:
	perl $0 fq1 <prefix> <Q cutoff for trim> <Max trim length> <#N allowed> <Q cutoff for filter> <min length of reads allowed> <Qual system(33|64)>
e.g	perl $0 sample.fq clean 30 10 1 30 35
USAGE
	exit;
}
sub openMethod {shift;return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")}

### BODY ###
my ($fq,$pfx,$Qt,$l,$n,$Qf,$lf,$Qsys,$debug) = @ARGV;
$Qsys ||= 33;

open FQ,"gzip -dc $fq |",or die "error\n";
open OUT,"|gzip >$pfx.clean.fq.gz" or die "error OUT\n";
open STAT,"> $pfx.clean.stat_out",or die "error\n";

my ($total, $remainQ, $sum_bp, $max_bp, $min_bp) = (0, 0, 0, 0, 1e9);
my %DEBUG if $debug;

while(<FQ>){
	#FQ info
	my ($seq,$num,$qual,$originLength,$Tlength,$len,$count,$avgQ) =();
	chomp;
	my @a =split;
	(my $fqID = $a[0]) =~ s/\/[12]$//;
	chomp($seq = <FQ>);
	chomp($num = <FQ>);
	chomp($qual= <FQ>);
	$total ++;
	$originLength = length($seq);
	# trim
	$Tlength = &Qstat($qual,$Qt,"trim",$l);
	$len = $originLength - $Tlength;
	$seq = substr($seq,0,$len);
	$qual= substr($qual,0,$len);
	# filter
	$count = $seq=~tr/N/N/;
	$avgQ  = &Qstat($qual,$Qf,"filter");
	if($count <= $n && $len >= $lf && $avgQ >= $Qf){		# N number & length limit judgement
		print OUT "$fqID length=$len\n$seq\n$num\n$qual\n";
		# stat
		$remainQ ++;
		$max_bp = ($max_bp > $len)?$max_bp:$len;
		$min_bp = ($min_bp < $len)?$min_bp:$len;
		$sum_bp += $len;
	} elsif($debug){
		#discard for N|discard for length|discard for PQ
		my $point =0;
		$point += ($count <= $n )?0:1;
		$point += ($len   >= $lf)?0:2;
		$point += ($avgQ  >= $Qf)?0:4;
		$DEBUG{$point} += 1;
	}
}
close FQ;
close OUT;

my $avgL = $sum_bp / $remainQ;
my $rate = $remainQ / $total;
my $tag = basename($pfx);
my $debugHead = ($debug)?"\tN>$n|Len<$lf|PQ<$Qf|N+Len|N+PQ|Len+PQ|HOMER":"";
for(1..7){$DEBUG{$_}||=0};
my $debugRes  = ($debug)?"\t$DEBUG{1}|$DEBUG{2}|$DEBUG{4}|$DEBUG{3}|$DEBUG{5}|$DEBUG{6}|$DEBUG{7}":"";

print STAT "Total\tmax\tmin\tavg\t#remain\trate\tSampleTAG(trim=$l,Qt=$Qt,N=$n,Qf=$Qf,min=$lf)$debugHead\n";
print STAT "$total\t$max_bp\t$min_bp\t$avgL\t$remainQ\t$rate\t$tag$debugRes\n";

close STAT;

# sub
sub Qstat {
	my ($q_l,$q_n,$method,$c_n) = @_;
	my $c = 0;

	if($method eq "trim"){
		for(my $i=length($q_l)-1;$i>=0;$i--){
			my $q=substr($q_l,$i,1);
			$q=ord$q;
			$q=$q - $Qsys;
			last if($q>=$q_n || $c>=$c_n);
			$c++;
		}
	}elsif($method eq "filter"){
		for(my $i=length($q_l)-1;$i>=0;$i--){
			my $q=substr($q_l,$i,1);
			$q=ord$q;
			$q=$q -$Qsys;
			$c += $q;
		}
		$c = $c/length($q_l);
		#$c = ($c/length($q_l) >= $q_n)?length($q_l):0;
	}
	return($c);
}



