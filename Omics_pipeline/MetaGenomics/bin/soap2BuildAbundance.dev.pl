#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Text;
use FindBin qw/$Bin/;

my ($sample_pair_1,$sample_pair_2,$sample_single,$dbs,$parameter,$ab,$ins,$sample,$prefix,$help,$workpath);
GetOptions(
        #"s:s"=>\$scaf_fa,
	#"g:s"=>\$gene_prediction_fa,
	"i1:s"=>\$sample_pair_1,
	"i2:s"=>\$sample_pair_2,
	"i3:s"=>\$sample_single,
	"DB:s" =>\$dbs,
	"par:s"=>\$parameter,
	"ab:i" =>\$ab,
	"ins:s"=>\$ins,
	"s:s"=>\$sample,
	"p:s"=>\$prefix,
	"o:s"=>\$workpath,
	"h:s"=>\$help,
        );

do{&usage;exit(1);} if ($help || !defined($prefix ||$sample_pair_1) );

$parameter ||= "m=226,x=426,r=0,l=30,M=4,S,p=8,v=5,S,c=0.95";
$parameter =~ s/^/ -/g; 
$parameter =~ s/,/ -/g;
$parameter =~ s/=/ /g;

$ab ||= 123;
$sample ||= $prefix;
$prefix ||= $sample;
chomp (my $pwd=`pwd`);
$workpath||=$pwd;

####SOAP####
my $soap_path = "$FindBin::Bin/soap2.22";
my $db_index1 = "/ifs1/ST_MD/USER/caixianghang/backup/MetaHit/27.1267sample_profile/list/db/4Group_uniqGene.div_1.fa.index";
my $db_index2 = "/ifs1/ST_MD/USER/caixianghang/backup/MetaHit/27.1267sample_profile/list/db/4Group_uniqGene.div_2.fa.index";
$dbs ||= join(',',$db_index1,$db_index2);
$dbs =~ s/,/ -D /g;
####conf####
`mkdir -p $workpath/$sample.gene.build`;

if($sample_single){
	my $cmd = "$soap_path -a $sample_single -D $dbs -o $workpath/$sample.gene.build/$prefix.soap.single.se $parameter 2>$workpath/$sample.gene.build/$prefix.soap.single.log\n";
	 print "$cmd";
	unless(system($cmd)){
		print "single soap finished!\n";
	}else{
		print "$?:single soap failed!\n";
		exit(($?>128?$?-128:$?));
	}
}

if($sample_pair_2){
	my $cmd = "$soap_path -a $sample_pair_1 -b $sample_pair_2 -D $dbs -o $workpath/$sample.gene.build/$prefix.soap.pair.pe -2 $workpath/$sample.gene.build/$prefix.soap.pair.se $parameter 2> $workpath/$sample.gene.build/$prefix.soap.pair.log\n";
	print "$cmd";
	unless(system($cmd)){
		print "pair soap finished!\n";
	}else{
		print "$?: pair soap failed!\n";
		exit(($?>128?$?-128:$?));
	}
}elsif($sample_pair_1){
	my $cmd = "$soap_path -a $sample_pair_1 -D $dbs -o $workpath/$sample.gene.build/$prefix.soap.SE.se $parameter 2> $workpath/$sample.gene.build/$prefix.soap.SE.log\n";
	print "$cmd";
	unless(system($cmd)){
		print "SE soap finished!\n";
	}else{
		die "ERR $?: SE soap failed!\n";
		#exit(($?>128?$?-128:$?));
	}
}

unless(system("gzip -f $workpath/$sample.gene.build/*.[ps]e")){
	exit(0);
}else{
	exit(($?>128?$?-128:$?));
}

########################
# sub function
########################
sub usage {
        print <<EOD;
Description: This program is used to produce IGC gene set profile.

Version 1.00 Feb 11,2014

Usage: perl $0 -i1 <reads_pair1.fq> -p <output file prefix> -o <output dir Name>

        Options:
	    -i1 <str>   the sample pair 1
	    -i2 <str>   the sample pair 2
	    -i3 <str>   the sample single 
		-par <str>	the parameters for soap, default as ",r=0,l=30,M=4,S,p=8,v=5,S,c=0.95" (if change, retain the first comma)
	    -ins <str>  the insrert size file, which contains sample name and insert size per row. e.g: Sample1	350
            -o  <str>   output dir Name
            -p  <str>   output file prefix
		-ab <int>	the abundance type u wanna build,default as "123" :
						1	reads abundance
						2	base abundance
						3	adjust abundance
	    -h  <str>   pod this help prefix
 data: 2015-3-26
 modified by linyuxiang\@genomics.cn
 Modified:	2016-01-07	fangchao\@genomics.cn

EOD
exit(1);
}

