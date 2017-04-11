#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Text;
use FindBin qw/$Bin/;

my ($sample_pair_1,$sample_pair_2,$sample_single,$ins,$prefix,$help,$workpath);
GetOptions(
        #"s:s"=>\$scaf_fa,
	#"g:s"=>\$gene_prediction_fa,
	"i1:s"=>\$sample_pair_1,
	"i2:s"=>\$sample_pair_2,
	"i3:s"=>\$sample_single,
	"ins:s"=>\$ins,
	"p:s"=>\$prefix,
	"o:s"=>\$workpath,
	"h:s"=>\$help,
        );

do{&usage;exit(1);} if ($help || !defined($prefix ||$sample_pair_1) );

chomp (my $pwd=`pwd`);
$workpath||=$pwd;

####SOAP####
my $soap_path = "/ifs1/ST_MD/USER/chenwn/bin/profiling/bin/soap2.22";
my $db_index1 = "/ifs1/ST_MD/USER/caixianghang/backup/MetaHit/27.1267sample_profile/list/db/4Group_uniqGene.div_1.fa.index";
my $db_index2 = "/ifs1/ST_MD/USER/caixianghang/backup/MetaHit/27.1267sample_profile/list/db/4Group_uniqGene.div_2.fa.index";
####conf####
`mkdir -p $workpath/$prefix.gene.build`;
if($sample_single){
	my $cmd = "$soap_path -a $sample_pair_1 -b $sample_pair_2 -D $db_index1 -D $db_index2 -m 226 -x 426 -o $workpath/$prefix.gene.build/$prefix.soap.pair.pe -2 $workpath/$prefix.gene.build/$prefix.soap.pair.se -r 0 -l 30 -M 4 -S -p 8 -v 5 -S -c 0.95 2> $workpath/$prefix.gene.build/$prefix.soap.pair.log\n";
	print STDERR "$cmd\n";
	`$cmd`;
	$cmd = "$soap_path -a $sample_single -D $db_index1 -D $db_index2 -o $workpath/$prefix.gene.build/$prefix.soap.single.se -r 0 -l 30 -M 4 -S -p 8 -v 5 -S -c 0.95 2>$workpath/$prefix.gene.build/$prefix.soap.single.log";
	 print STDERR "$cmd\n";
	`$cmd`;
}
elsif($sample_pair_2){
	my $cmd = "$soap_path -a $sample_pair_1 -b $sample_pair_2 -D $db_index1 -D $db_index2 -m 226 -x 426 -o $workpath/$prefix.gene.build/$prefix.soap.pair.pe -2 $workpath/$prefix.gene.build/$prefix.soap.pair.se -r 0 -l 30 -M 4 -S -p 8 -v 5 -S -c 0.95 2> $workpath/$prefix.gene.build/$prefix.soap.pair.log";
	print STDERR "$cmd\n";
	`$cmd`;
	}
elsif($sample_pair_1){
	my $cmd = "$soap_path -a $sample_pair_1 -D $db_index1 -D $db_index2 -m 226 -x 426 -o $workpath/$prefix.gene.build/$prefix.soap.pair.pe -2 $workpath/$prefix.gene.build/$prefix.soap.pair.se -r 0 -l 30 -M 4 -S -p 8 -v 5 -S -c 0.95 2> $workpath/$prefix.gene.build/$prefix.soap.pair.log\n";
	print STDERR "$cmd";
	`$cmd`;print STDERR "soap finished!\n";
}
chomp(my $f = `ls $workpath/$prefix.gene.build/*.[ps]e|head -1`);
chomp(my $g = `ls $workpath/$prefix.gene.build/*.[ps]e.gz|head -1`);
if (-e $f){
	`find $workpath/$prefix.gene.build/*.[ps]e >$workpath/$prefix.soap.list`;
}elsif(-e $g){
	`find $workpath/$prefix.gene.build/*gz > $workpath/$prefix.soap.list`;
}else{
	die "can not find soap result files.".$!;
}


my $cmd = "perl /ifs1/ST_MD/USER/chenwn/bin/profiling/bin/gene_Profiling.pl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len $ins $workpath/$prefix.soap.list $workpath/$prefix\n";
`$cmd`;print STDERR "reads abudance file built\n";
$cmd = "perl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/profile/05.LostReads_ProfileAdjust/gene_Profiling.cxh.v3.all.pl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len $ins $workpath/$prefix.soap.list $workpath/$prefix.base\n";
#print STDERR $cmd;`$cmd`;
`$cmd`;print STDERR "base abudance file built\n";
$cmd = "perl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/profile/05.LostReads_ProfileAdjust/gene_Profiling.cxh.v4.all.pl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len $ins $workpath/$prefix.soap.list $workpath/$prefix.aj\n";
#print STDERR $cmd;`$cmd`;
`$cmd`;print STDERR "adjusted abudance file built\n";

####### compress the results to save space ### added by fangchao@genomics.cn
`gzip $workpath/$prefix.gene.build/*.[ps]e`;
`gzip -f $workpath/$prefix*.abundance`;

sub usage {
        print <<EOD;
Description: This program is used to produce IGC gene set profile.

Version 1.00 Feb 11,2014

Usage: perl $0 -i1 <reads_pair1.fq> -p <output file prefix> -o <output dir Name>

        Options:
	    -i1 <str>   the sample pair 1
	    -i2 <str>   the sample pair 2
	    -i3 <str>   the sample single 
	    -ins <str>  the insrert size file, which contains sample name and insert size per row. e.g: Sample1	350
            -o  <str>   output dir Name
            -p  <str>   output file prefix
	    -h  <str>   pod this help prefix
 data: 2015-3-26
 modified by linyuxiang\@genomics.cn
EOD
exit(1);
}

