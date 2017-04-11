#!/usr/bin/perl -w
use strict;
################################################################################
unless(4==@ARGV) {
    &usage;
    exit;
}
################################################################################
my($list_f,$order_f,$row,$out) = @ARGV;
my(@order,%list,@info,$i,%class,%cover,%Count,%sum_abun);
################################################################################
open IN,"<$list_f" || die "read $list_f $!\n";
while(<IN>) {
    chomp;
    @info=split /\//;
    $info[-1]=~/^(\S+)\.(abun\.GC\.readsGC|readsGC)(|\.gz)$/;
	$info[-1]=~/^([^\.]+)(\S+)$/;
    my $name =$1;
    $list{$name}=$_;
}
close IN;
################################################################################
open IN,$order_f or die "read $order_f $!\n";
while(<IN>) {
    chomp;
    push(@order,$_);
}
close IN;
################################################################################
for($i=0;$i<@order;++$i) {
	my $openMethod = ($list{$order[$i]} =~ /\.gz$/)?"gzip -dc $list{$order[$i]} |":"$list{$order[$i]}";
	open IN,$openMethod  or die "$!\n";
#    <IN>;		# MAKE SURE YOUR ABUNDANCE FILE GOT A HEADER!!!
    while(<IN>) {
        chomp;
        @info=split /\t/;
		if ($.==1) {
			next unless $info[$row] =~ /\d/;
		}
        $class{$info[0]}.="\t".$info[$row];
		$cover{$info[0]} ||= 0; 
		$Count{$i} ||= 0;
		if ($info[$row] > 0){ $cover{$info[0]} ++ ; $Count{$i} ++ ;$sum_abun{$info[0]} += $info[$row] };
    }
    close IN;
}
################################################################################
open OT,">$out.profile" or die "write $out $!\n";
open CrT,">$out.cover" || die $!;
open CtT,">$out.count" || die $!;
for($i=0;$i<@order;++$i) {
    print OT "\t",$order[$i];
	print CtT "$order[$i]\t$Count{$i}\n";
}
print OT "\n";
foreach $i(sort {$a<=>$b} keys %class) {
	my $avg;
	if($cover{$i} >0){
		$avg = $sum_abun{$i} / $cover{$i};
	}else{$avg = 0}
    print OT $i,$class{$i},"\n";
	print CrT "$i\t$cover{$i}\t$avg\n";
}
close OT;
close CrT;
################################################################################
sub usage {
    print STDERR<<USAGE;
Description
	This programme is to combine profiling table
Usage:
	perl $0 [file.list] [order] [row] [outfile prefix]
		Row must be in 1(pairs),2(base_abundance),3(reads_abundance),4(depth_abundance)

	Author Libranjie,zhouyuanjie\@genomics.org.cn
	updated 2014/12/5 linyuxiang\@genomics.cn
	customized by fangchao\@genomics.cn # 20160118
USAGE
}
