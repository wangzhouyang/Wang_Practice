#!/usr/bin/perl -w
use warnings;
use strict;
use File::Basename; 

die &usage if @ARGV <=7;
sub usage {
	print <<USAGE;
usage:
pe pattern:
#	perl $0 fq1,fq2 <outdir> <qt> <limit> <N num> <qf> <lf> <PhQ> <discard>
e.g	perl $0 sample_1.fq,sample_2.fq clean 20 10 1 15 0 
		dir	    path/to/output/
		qt		PhredQual cutoff for trim from tail
		limit	maximum trim length
		N num	N base tolerance
		qf		PhredQual cutoff for filter. The average Qual lower than qf will be discarded.
		lf		minimun length of trimed reads.
		PhQ     Qual system(33|64)
        discard output discard reads  
USAGE
	exit;
}
sub openMethod {$_=shift;return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")}

my ($fq,$out,$qt,$l,$n,$qf,$lf,$PhQ,$discard,$suffix) = @ARGV;
my @fqs = split(",",$fq);
my ($fq1,$fq2) = @fqs;
open FQ1, &openMethod($fq1) or die "error\n";
open FQ2, &openMethod($fq2) or die "error\n";
##
$suffix ||= 'fq.gz';
my $method=($suffix eq "fq")?'>':'|gzip >';

open OUT1,"$method $out.clean.1.$suffix" or die "error OUT1\n";
open OUT2,"$method $out.clean.2.$suffix" or die "error OUT2\n";
open OUT3,"$method $out.clean.single.$suffix" or die "error OUT3\n";
if($discard){
    open DIS,"|gzip >$out.clean.dirty.fq.gz" or die "error DIS\n";
}
open STAT,"> $out.clean.stat_out",or die "error\n";

my @total = (0,0);
my(@remainN,@remainQ,@sum_bp)= ();
my @max_bp = (0,0);
my @min_bp = (10000,10000);

my %READ;
my($single1,$single2,$paired) = (0,0,0);
while(<FQ1>){
    my $flag = 0;
    my ($readsSEQ1,$readsSEQ2)=(0,0);
    #FQ1
    my ($seq,$num,$quality,$originLength,$Tscore,$len,$count,$Qscore) =();
    chomp;
    my @a =split;
    (my $out1 = $a[0]) =~ s/\/[12]$//;
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
    $Qscore = &Qstat($quality,$qf,"filter");
    $readsSEQ1 = "$out1/1 length=$len\n$seq\n$num\n$quality\n";
    if($count <= $n && $len >= $lf && $Qscore >= $qf){		# N number & length limit judgement
        $flag += 1;
        $remainN[0] ++ ; 
        # filter more
        $remainQ[0] ++;
        $max_bp[0] = ($max_bp[0] > $len)?$max_bp[0]:$len;
        $min_bp[0] = ($min_bp[0] < $len)?$min_bp[0]:$len;
        $sum_bp[0] += $len;
    }
    # FQ2
    ($seq,$num,$quality,$originLength,$Tscore,$len,$count,$Qscore) =();
    chomp($fq2 = <FQ2>);
    my @b= split(/\s|\t/,$fq2);
    (my $out2 = $b[0]) =~ s/\/[12]$//;
    print STDERR "$out1/1 \& $out2/2 fq id isn't paired ." if ( $out1 ne $out2);
    chomp($seq=<FQ2>);
    chomp($num=<FQ2>);
    chomp($quality=<FQ2>);
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
    $readsSEQ2 = "$out2\/2 length=$len\n$seq\n$num\n$quality\n";
    if($count <= $n && $len >= $lf && $Qscore >= $qf){
        $flag += 2;
        $remainN[1] ++ ;
        $remainQ[1] ++;
        $max_bp[1] = ($max_bp[1] > $len)?$max_bp[1]:$len;
        $min_bp[1] = ($min_bp[1] < $len)?$min_bp[1]:$len;
        $sum_bp[1] += $len;
    }

    if($flag == 1){
        $single1 += 1;
        print OUT3 $readsSEQ1;
        print DIS $readsSEQ2 if $discard;
    }elsif($flag == 2){
        $single2 += 1;
        print DIS $readsSEQ1 if $discard;
        print OUT3 $readsSEQ2;
    }elsif($flag == 3){
        $paired += 1;
        print OUT1 $readsSEQ1;
        print OUT2 $readsSEQ2;
    }
}
close FQ1;
close FQ2;
close OUT1;
close OUT2;
close OUT3;
close DIS if $discard;

my $avg1 = $sum_bp[0] / $total[0];
my $avg2 = $sum_bp[1] / $total[1] if $total[1] >0 ;
my $rate1 = $remainQ[0] / $total[0];
my $rate2 = $remainQ[1] / $total[1] if $total[1] >0 ;
my $tag = basename($out);

print STAT "Total\tmax1\tmin1\tmax2\tmin2\tremain_pair\tremain1\tremain2\tSampleTAG(trim_limit=$l,Qt=$qt,N=$n,Qf=$qf,min=$lf)\n";
print STAT "$total[0]\t$max_bp[0]\t$min_bp[0]\t$max_bp[1]\t$min_bp[1]\t$paired\t$single1\t$single2\t$tag\n";

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
            $q=$q - $PhQ;
            $c += $q
        }
        $c = $c / length($q_l);
    }elsif($method eq "trim"){
        for(my $i=length($q_l)-1;$i>=0;$i--){
            my $q=substr($q_l,$i,1);
            $q=ord$q;
            $q=$q - $PhQ;
            last if($q>=$q_n || $c>=$c_n);
            $c++;
        }
    }
    return($c);
}



