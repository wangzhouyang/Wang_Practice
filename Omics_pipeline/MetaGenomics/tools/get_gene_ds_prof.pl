#!/usr/bin/perl -w
use strict;
unless(@ARGV >= 4) {
	&usage;
	exit;
}

my($reflist,$insertsize,$soaplistlist,$ds,$out)=@ARGV;

my(%lengths,%insert,@temp,$i,$j,%abundance,$totalreads,$totalabundance,@genename,%ids);
###########################################################
############### reading referencr length
if($reflist =~ /gz$/){
    open IN,"gzip -dc $reflist |" or die "$reflist $!\n";
}else{
    open IN,$reflist or die "$reflist $!\n";
}
while(<IN>) {
	chomp;
	@temp=split;
	next if ($temp[2] == 0);
    $lengths{$temp[1]}=$temp[2];
	$ids{$temp[1]}=$temp[0];
	push (@genename,$temp[1]);
}
close IN;

################ reading insert size
if ($insertsize ne "SE"){
	open IN,$insertsize or die "$insertsize $!\n";
	while(<IN>) {
		chomp;
		@temp=split;
		$insert{$temp[0]}=$temp[1];
	}
}

################ reading list of soap's result
open IN,$soaplistlist or die "$soaplistlist $!\n";
my %readsnum;
while(<IN>) {
	chomp;
	################ each is soap result list
	@temp=split /\//,$_;
	$temp[-1] =~ /(\S+)\.soap.*/;
	$j=$1;
#	if($j=~/^(\S+\.\D+)(\d+)\-\d(\_\d+)/) {
#		$j=$1."\-".$2.$3;
#	}
	my $tmp=$_;
    if ($insertsize ne "SE" && not exists $insert{$j}) {
        die "Unknown InsertSize of sample $j,\nplease check!\n";
    }
	&stat_readsnum($tmp,$insert{$j});
	last if $ds >0 && $totalreads >= $ds;
}
close IN;

################ print readsnum and calculate abundance

for($i=0;$i<@genename;++$i) {
	unless(defined $readsnum{$genename[$i]}) {
		$readsnum{$genename[$i]}=0;
	}
	$abundance{$genename[$i]}=$readsnum{$genename[$i]}/$lengths{$genename[$i]};
	$totalabundance+=$abundance{$genename[$i]};
	#$totalreads+=$readsnum{$genename[$i]};
}
#print "total_reads:\t$totalreads\n";
#print "total_abundance:\t$totalabundance\n";

################ print abundance 
open OT,">$out.abundance" or die "$out.abundance $!\n";
open ST,">$out.abundance.size" or die "cannot wirte to $out.abundance.size.\n $!\n";
print ST "$totalreads\n";close ST;
print OT "ID\treads_pairs\tbase_abundance\treads_abundance\n";
for($i=0;$i<@genename;++$i) {
	print OT
		"$ids{$genename[$i]}\t$readsnum{$genename[$i]}\t",$readsnum{$genename[$i]}/$totalreads,"\t",$abundance{$genename[$i]}/$totalabundance,"\n";
}
close OT;

###########################################################
################ submit stat reads number function
sub stat_readsnum{
	my($in,$inser)=@_;
	if($in =~ /\.pe$/ or $in =~ /\.pe\.gz$/ or $insertsize eq "SE") {
		if($in =~ /\.[ps]e$/) {
			open PE,$in or die "$in $!\n";
		}elsif($in =~ /\.[ps]e\.gz$/) {
			open PE,"gzip -dc $in|" or die "$in $!\n";
		}
		my @tmp;
		my $count_per_read = ($insertsize eq "SE")?"1":"0.5";
		while(<PE>) {
			chomp;
			@tmp=split;
			if($tmp[3] ==1) {
				$readsnum{$tmp[7]}+= $count_per_read;
				$totalreads += $count_per_read;
				last if $ds >0 && $totalreads >= $ds;
			}
		}
		close PE;
	}elsif($in =~ /\.se$/ or $in =~ /\.se\.gz$/) {
		if($in =~ /\.se$/) { 
			open SE,$in or die "$in $!\n";
		}elsif($in =~ /\.se\.gz$/) {
			open SE,"gzip -dc $in|" or die "$in $!\n";
		}
		my(@tmp,%past);
		while(<SE>) {
			chomp;
			@tmp=split;
			if($tmp[3] ==1) {
				if($tmp[6] eq '+') {
					if(($lengths{$tmp[7]}-$tmp[8])<$inser+100) {
						$tmp[0]=~/^(\S+)\/[12]/;
						$past{$tmp[7]}{$1}=1;
					}
				} elsif($tmp[6] eq '-') {
					if(($tmp[8])<$inser-$tmp[5]+100) {
						$tmp[0]=~/^(\S+)\/[12]/;
						$past{$tmp[7]}{$1}=1;
					}
				}
			}
		}
		close SE;
		foreach $i(sort keys %past) {
			foreach my $read(sort keys %{$past{$i}}) {
				if($past{$i}{$read}) {
					++$readsnum{$i};
				}
			}
		}
	}
}

################ submit usage function
sub usage{
	print STDERR "Description
        This programme is to creat profiling result!
        Usage:$0 [reference.length.list] [insertsize.list] [soap.list] [outfile-prefix]
        In this programme
		[insertsize.list]	a couple of paired reads in pair end is one pair,
					a couple of single end only within insert size is one pair.
					But SINGLE-END reads DO NOT exist a insert issue, so just enter \"SE\" when you perform a single-end reads abundance.
        5th verison by Wed Jul 21 09:27:02 2010
        Author Libranjie,zhouyuanjie\@genomics.org.cn
		Modifier fangchao\@genomics.cn\n"
}
