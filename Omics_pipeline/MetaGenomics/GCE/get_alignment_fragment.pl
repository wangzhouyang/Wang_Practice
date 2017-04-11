#!/usr/bin/perl -w
#use strict;
unless(3==@ARGV) {
	&usage;
	exit;
}

my($reflist,$soaplistlist,$out)=@ARGV;

my($i,$j,%ids,%pos,%open);
###########################################################
sub openMethod {my $f=shift;my $m=($f=~/gz$/)?"gzip -dc $f|":"$f";return($m);}
system("mkdir -p $out");
############### reading referencr length
if($reflist =~ /gz$/){
    open IN,"gzip -dc $reflist |" or die "$reflist $!\n";
}else{
    open IN,$reflist or die "$reflist $!\n";
}
while(<IN>) {
	chomp;
	my @temp=split;
	$temp[4] =~ s/\s/\_/g;
	$temp[4] =~ s/\//\_OR\_/g;
	$ids{$temp[3]}=$temp[4];		#$ids{GENE NAME}= SPECIES
	$pos{$temp[3]}=$temp[2];
}
close IN;

################ reading list of soap's result
open IN,$soaplistlist or die "$soaplistlist $!\n";
my %readsnum;
while(<IN>) {
	chomp;
	################ each is soap result list
	my @temp=split /\//,$_;
	$temp[-1] =~ /(\S+)\.soap.*/;
	$j=$1;
#	if($j=~/^(\S+\.\D+)(\d+)\-\d(\_\d+)/) {
#		$j=$1."\-".$2.$3;
#	}
	my $tmp=$_;
	
	open PE,&openMethod($tmp) or die "$tmp $!\n";
	while(<PE>) {
		chomp;
		my @temp=split;
		if($temp[3] ==1) {
			next if not defined $pos{$temp[7]};
			my $locate = $temp[8] + $pos{$temp[7]};
			my $sp = $ids{$temp[7]};
			my $op = "POS".$sp;
			if (not defined $open{$op}){
				$open{$op} = $op;
				open $op,"> $out/$sp.pos" or die "cant open $out/$sp.pos, $!\n";
			}
			print $op "$sp\t$locate\t$temp[5]\n";
		}
	}
	close PE;
}

foreach my $op (values %open){close "POS".$op;}

################ submit usage function
sub usage{
	print STDERR "Description
        Usage:$0 [reference.length.list] [soap.list] [outfile-prefix]
		Modifier fangchao\@genomics.cn\n"
}
