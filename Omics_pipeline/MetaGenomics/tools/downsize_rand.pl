#!/usr/bin/perl -w
use strict;
die &usage unless @ARGV >= 4;
my($originSize,$targetSize,$reads,$length,$out,$single_double) = @ARGV;
# $single_double means weather the reads counts need to be counted tiwce (because some pair-end abundance calculate 1 read as half a *pair*)
$single_double ||= 1;
sub openMethod{my $f=shift;return((($f =~ /\.gz$/)?"gzip -dc $f|":"$f"))}
sub writeM{my $f=shift;return((($f =~ /\.gz$/)?"|gzip > $f":"> $f"))}
################################## downsize ######################################
srand($originSize);
my (%reads_num,@rands);
open I,&openMethod($reads) or die "$!\n";
open ST, &writeM($out) or die "$!\n";
my $total_reads=$originSize;
my $k=0;
while(<I>)
{
    chomp;
    my @temp=split;
	next if $.==1 && $temp[1]=~/[a-z]/;
	$total_reads+= $temp[1] * $single_double if $originSize==0;
    $reads_num{$temp[0]}=($originSize-2*$targetSize<0)?$temp[1] * $single_double:0;
	if ($temp[1] >0){
		my $t=$temp[1] * $single_double;
		while($t){
			$rands[$k] = $temp[0];
			$k ++;$t --;
		}
	}
}
close I;

my $time=$total_reads-$targetSize; 
if ($time < 0){
	print ST "Less than $targetSize ($total_reads). Discarded.\n";
	exit;
}else{
	$time = ($time < $targetSize)?$time:$targetSize;
	my $one  = ($time < $targetSize)?-1:1;
	while($time){
	  	my $point=int(rand(@rands));
	    $reads_num{$rands[$point]} += $one;
		splice(@rands,$point,1);
	    $time --;
	}
}
	

################################ abundance ################################

open I,&openMethod($length) or die "$!\n";
my $gene_n=0;###size of gene length profile
my %reads_abundance;### reads_abundance=reads number/gene length
my $total_abundance=0;### total_abundance of reads_abundance
### base_abundance=reads_num/total_reads_number 
while(<I>){  
    chomp;
    my @temp=split "\t",$_;
    $gene_n++;
	$reads_num{$temp[0]}||=0;
	if($reads_num{$temp[0]}==0){
		$reads_abundance{$temp[0]}=0;
	}else{
		if($temp[2]==0){next;}else{
			$reads_abundance{$temp[0]}=$reads_num{$temp[0]}/$temp[2];
			$total_abundance+=$reads_abundance{$temp[0]};
		}	
	}
}
close I;
print ST "ID\treads_pairs\tbase_abundance\treads_abundance\n";
for my $i(1 .. $gene_n) {
	print ST "$i\t$reads_num{$i}\t",$reads_num{$i}/$targetSize,"\t",$reads_abundance{$i}/$total_abundance,"\n";
}

close ST;
sub usage
{
    print "usage:perl $0 [origin read number] [targetSize] [input file] [len.info file] [output ifle]\n
			The number of reads input contain. [0] if unknown
			The number of reads you wanna trim.
			A list contain the gene id and reads number(gene & reads number before downsize)
			A list contain the gene id and its length(gene length reference)
			output prefix\n";
exit;
}
