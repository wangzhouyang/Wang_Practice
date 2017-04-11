#!/usr/bin/perl -w
use warnings;
use strict;
use File::Basename; 

die &usage if @ARGV <=7;
sub usage {
	print <<USAGE;
usage:
pe pattern:
#	perl $0 fq1,fq2 <prefix> <qt> <limit> <N num> <qf> <lf>
se pattern:
	perl $0 fq1 <prefix> <qt> <limit> <N num> <qf> <lf> <avgQ>(Average Qual cutoff mode, option)
e.g	perl $0 sample.fq clean 20 10 1 15 0 
		prefix	path/to/output/file_prefix
		qt		PhredQual cutoff for trim from tail
		limit	maximum trim length
		N num	N base tolerance
		qf		PhredQual cutoff for filter. The average Qual lower than qf will be discarded.
		lf		minimun length of trimed reads.
		avgQ	
USAGE
	exit;
}

my ($fq,$pfx,$qt,$l,$n,$qf,$lf,$avgQ) = @ARGV;
my @fqs = split(",",$fq);
my ($fq1,$fq2) = (@fqs == 2)?@fqs:($fq,$fq);
if($fq1 eq $fq2){
	open FQ1,"gzip -dc $fq1 |",or die "error\n";
	open OUT3,"|gzip >$pfx.clean.fq.gz",or die "error\n";
}else{
	open FQ1,"gzip -dc $fq1 |",or die "error\n";
	open FQ2,"gzip -dc $fq2 |",or die "error\n";
	open OUT1,"|gzip >$pfx.clean.1.fq.gz" or die "error OUT1\n";
	open OUT2,"|gzip >$pfx.clean.2.fq.gz" or die "error OUT2\n";
	open OUT3,"|gzip >$pfx.clean.single.fq.gz" or die "error OUT3\n";
}
open STAT,"> $pfx.clean.stat_out",or die "error\n";

my @total = (0,0);
my(@remainN,@remainQ,@sum_bp)= ();
my @max_bp = (0,0);
my @min_bp = (10000,10000);

my %READ;
while(<FQ1>){
	my $flag = 0;
	my ($readsSEQ1,$readsSEQ2)=(0,0);
	#FQ1
	my ($seq,$num,$quality,$originLength,$Tscore,$len,$count,$Qscore) =();
	chomp;
	my @a =split;
	(my $pfx1 = $a[0]) =~ s/\/[12]$//;
	chomp($seq=<FQ1>);
	chomp($num=<FQ1>);
	chomp($quality=<FQ1>);
	$total[0] ++;
	$originLength = length($seq);
	# trim
	$Tscore = &Qstat($quality,$qt,"trim",$l);
	$len = $originLength-$Tscore;
	$seq=substr($seq,0,$len);
	$quality=substr($quality,0,$len);
	# filter
	$count=$seq=~tr/N/N/;
	if($avgQ){
		$Qscore = &Qstat($quality,$qf,"avgQ");
	}else{
		$Qscore = &Qstat($quality,$qf,"filter");
	}
	if($count <= $n && $len >= $lf && ($Qscore*2) <= length($quality)){		# N number & length limit judgement
		$flag += 1;
		$remainN[0] ++ ; 
		# filter more
		my $PFX = ( $fq1 eq $fq2)?$pfx1:"$pfx1\/1";
		$readsSEQ1 = "$PFX length=$len\n$seq\n$num\n$quality\n";
		#print OUT1 $readsSEQ;
		$remainQ[0] ++;
		$max_bp[0] = ($max_bp[0] > $len)?$max_bp[0]:$len;
		$min_bp[0] = ($min_bp[0] < $len)?$min_bp[0]:$len;
		$sum_bp[0] += $len;
	}
	# FQ2
	($seq,$num,$quality,$originLength,$Tscore,$len,$count,$Qscore) =();
	if( $fq1 ne $fq2){

		chomp($fq2 = <FQ2>);
		my @b= split(/\s|\t/,$fq2);
		(my $pfx2 = $b[0]) =~ s/\/2$//;
		print STDERR "$pfx1/1 & $pfx2/2 fq id isn't paired ." if ( $pfx1 ne $pfx2);
		chomp(my $seq=<FQ2>);
		chomp(my $num=<FQ2>);
		chomp(my $quality=<FQ2>);
		$total[1] ++;
		$originLength = length($seq);
		# trim
		$Tscore = &Qstat($quality,$qt,"trim",$l);
		$len = $originLength-$Tscore;
		$seq=substr($seq,0,$len);
		$quality=substr($quality,0,$len);
		# filter
		$count=$seq=~tr/N/N/;
		$Qscore = &Qstat($quality,$qf,"filter");
		if($count <= $n && $len >= $lf && ($Qscore*2) <= length($quality)){
			$flag += 2;
			$remainN[1] ++ ;
			$readsSEQ2 = "$pfx2\/2 length=$len\n$seq\n$num\n$quality\n";
			$remainQ[1] ++;
			$max_bp[1] = ($max_bp[1] > $len)?$max_bp[1]:$len;
			$min_bp[1] = ($min_bp[1] < $len)?$min_bp[1]:$len;
			$sum_bp[1] += $len;
		}
	}
	# print 
	if($flag == 1){
		print OUT3 $readsSEQ1;
#		print OUT3 $readsSEQ2 if $fq1 ne $fq2;
	}elsif($flag == 2){
#		print OUT2 $readsSEQ2;
		print OUT3 $readsSEQ2;
	}elsif($flag == 3){
		print OUT1 $readsSEQ1;
		print OUT2 $readsSEQ2;
	}
}
close FQ1;
close FQ2;
close OUT1;
close OUT2;
close OUT3;

my $avg1 = $sum_bp[0] / $total[0];
my $avg2 = $sum_bp[1] / $total[1] if $total[1] >0 ;
my $rate1 = $remainQ[0] / $total[0];
my $rate2 = $remainQ[1] / $total[1] if $total[1] >0 ;
my $tag = basename($pfx);
unless ($fq1 eq $fq2){
	print STAT "Total\tmax1\tmin1\tavg1\tmax2\tmin2\tavg2\tremain1\tremain2\trate1\trate2\tSampleTAG(trim_limit=$l,Qt=$qt,N=$n,Qf=$qf,min=$lf)\n";
	print STAT "$total[0]\t$max_bp[0]\t$min_bp[0]\t$avg1\t$max_bp[1]\t$min_bp[1]\t$avg2\t$remainQ[0]\t$remainQ[1]\t$rate1\t$rate2\t$tag\n";
}else{
	print STAT "Total\tmax\tmin\tavg\tremain\trate\tSampleTAG(trim_limit=$l,Qt=$qt,N=$n,Qf=$qf,min=$lf)\n";
	print STAT "$total[0]\t$max_bp[0]\t$min_bp[0]\t$avg1\t$remainQ[0]\t$rate1\t$tag\n";
}

close STAT;

# sub
sub lengthAvg{
	my($l,$max,$min,$sum) = @_;
	$max = ($max > $l)?$max:$l;
	$min = ($min < $l)?$min:$l;
	$sum += $l;
	return($max,$min,$sum);
}

sub Qstat {
	my ($q_l,$q_n,$method,$c_n) = @_;
	my $c = 0;
	if ($method eq "filter"){
		for(my $i=0;$i<length($q_l);$i++){
			my $q=substr($q_l,$i,1);
			$q=ord$q;
			$q=$q-33;
			if($q<=$q_n){$c++};
		}
	}elsif($method eq "trim"){
		for(my $i=length($q_l)-1;$i>=0;$i--){
			my $q=substr($q_l,$i,1);
			$q=ord$q;
			$q=$q-33;
			last if($q>=$q_n || $c>=$c_n);
			$c++;
		}
	}elsif($method eq "avgQ"){
		for(my $i=length($q_l)-1;$i>=0;$i--){
			my $q=substr($q_l,$i,1);
			$q=ord$q;
			$q=$q-33;
			$c += $q;
		}
		$c = ($c/length($q_l) >= $q_n)?0:length($q_l);
	}
	return($c);
}



