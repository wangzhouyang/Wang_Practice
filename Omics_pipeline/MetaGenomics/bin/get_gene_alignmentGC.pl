#!/usr/bin/perl -w
use strict;
unless(5==@ARGV) {
	&usage;
	exit;
}

my($reflist,$insertsize,$soaplistlist,$GCfileList,$out)=@ARGV;

my(%lengths,%insert,@temp,$i,$j,%abundance,$totalreads,$totalabundance,@genename,%ids,%GClen,%readsGC,%readslen,%readsid);
###########################################################
sub openMethod {my $f=shift;my $m=($f=~/gz$/)?"gzip -dc $f|":"$f";return($m);}
############### reading referencr length
$reflist=($reflist eq "default")?"/ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len":$reflist;
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
	push (@genename,$temp[0]);
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
open GL,$GCfileList or die "$GCfileList $!\n";
my %readsnum;
while(<IN>) {
	chomp;
	my $GCline = <GL>; chomp($GCline);
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
	&hash_GC($GCline);
#	undef %GClen;
}
close GL;
close IN;

################ print readsnum and calculate abundance
=cut
for($i=0;$i<@genename;++$i) {
	unless(defined $readsnum{$genename[$i]}) {
		$readsnum{$genename[$i]}=0;
		$readsGC{$genename[$i]}{'GC'}=0;
		$readsGC{$genename[$i]}{'len'}=1;
	}
	$GCcontent{$genename[$i]}=$readsGC{$genename[$i]}{'GC'}/$readsGC{$genename[$i]}{'len'};
	$abundance{$genename[$i]}=$readsnum{$genename[$i]}/$lengths{$genename[$i]};
	$totalabundance+=$abundance{$genename[$i]};
	$totalreads+=$readsnum{$genename[$i]};
}
=cut
#print "total_reads:\t$totalreads\n";
#print "total_abundance:\t$totalabundance\n";

################ print abundance 
open OT,">$out.readsGC" or die "$out.readsGC $!\n";
open ST,">$out.readsGC.size" or die "cannot wirte to $out.readsGC.size.\n $!\n";
print OT "ID\treads_pairs\tGC_content\n";
for($i=0;$i<@genename;++$i) {
	unless(defined $readsnum{$genename[$i]}){
		$readsnum{$genename[$i]}=0;
		$readsGC{$genename[$i]}=0;
		$readslen{$genename[$i]}=1;
	}
	unless(defined $readsGC{$genename[$i]}){
#		$readsGC{$genename[$i]}=0;$readslen{$genename[$i]}=1;
		die "Gene $genename[$i] got $readsnum{$genename[$i]} reads but miss GC info!\n";
	}

	$totalreads += $readsnum{$genename[$i]};
	my $GCrate = $readsGC{$genename[$i]}/$readslen{$genename[$i]};
	print OT
#		"$ids{$genename[$i]}\t$readsnum{$genename[$i]}\t$GCrate\n";
		"$genename[$i]\t$readsnum{$genename[$i]}\t$GCrate\n";
}
close OT;
print ST "$totalreads\n";close ST; 
###########################################################
sub hash_GC{
	my $in = shift;
	open GC,&openMethod($in) or die "$in cannot open. $!\n";
	while(<GC>){
		chomp;
		my @a = split;
		next unless defined $readsid{$a[0]};
		$readslen{$readsid{$a[0]}} += $a[1];
		$readsGC{$readsid{$a[0]}} += $a[2];
#		$GClen{$a[0]} = "$a[1]\t$a[2]";
#		$GCGC{$a[0]} = $a[2];
	}
	close GC;
}
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
#				die "$tmp[0] miss GC info!\n" unless (defined $GClen{$tmp[0]});
#				unless (defined $GClen{$tmp[0]}) { $GClen{$tmp[0]} = "50\t25"};
#				my @gcs = split(/\t/,$GClen{$tmp[0]});

#				$readsnum{$tmp[7]} += $count_per_read;
				$readsnum{$ids{$tmp[7]}} +=$count_per_read;				
				$readsid{$tmp[0]} = $ids{$tmp[7]};
#				$readslen{$tmp[7]} += $gcs[0];
#				$readsGC{$tmp[7]} += $gcs[1];
#				delete $GClen{$tmp[0]};
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
        Usage:$0 [reference.length.list] [insertsize.list] [soap.list] [GC.list] [outfile-prefix]
        In this programme
		[insertsize.list]	a couple of paired reads in pair end is one pair,
					a couple of single end only within insert size is one pair.
					GC list should be in the same order with soap.list
					But SINGLE-END reads DO NOT exist a insert issue, so just enter \"SE\" when you perform a single-end reads abundance.
        5th verison by Wed Jul 21 09:27:02 2010
        Author Libranjie,zhouyuanjie\@genomics.org.cn
		Modifier fangchao\@genomics.cn\n"
}
