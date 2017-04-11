#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Text;
use FindBin qw/$Bin/;

my ($ab,$gl,$ins,$list,$prefix,$help,$workpath);
GetOptions(
	"ab:i" =>\$ab,
	"gl:s" =>\$gl,
	"ins:s"=>\$ins,
	"l:s"=>\$list,
	"p:s"=>\$prefix,
	"o:s"=>\$workpath,
	"h:s"=>\$help,
        );

$ab ||= 1;
chomp (my $pwd=`pwd`);
$workpath||=$pwd;
&usage unless ($gl && $ins && $list);

if ($ab =~ /1/){
#	my $cmd = "perl /ifs1/ST_MD/USER/chenwn/bin/profiling/bin/gene_Profiling.pl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len $ins $list $workpath/$prefix\n";
#	my $cmd = "perl $Bin/get_gene_profile.pl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len $ins $list $workpath/$prefix\n";
	my $cmd = "perl $Bin/get_gene_profile.pl $gl $ins $list $workpath/$prefix\n";
	print STDERR $cmd;`$cmd`;
	print STDERR "reads abundance file built\n";
}
if ($ab =~ /2/){
	my $cmd = "perl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/profile/05.LostReads_ProfileAdjust/gene_Profiling.cxh.v3.all.pl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len $ins $list $workpath/$prefix.base\n";
	print STDERR $cmd;`$cmd`;
	print STDERR "base abundance file built\n";
}
if ($ab =~ /3/){
	my $cmd = "perl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/profile/05.LostReads_ProfileAdjust/gene_Profiling.cxh.v4.all.pl /ifs1/ST_MD/PMO/SZC08004_MetaHIT/User/caixianghang/06.Profile/1.GeneProfile/list/760MetaHit_139HMP_368PKU_511Bac.uniq.fa.len $ins $list $workpath/$prefix.aj\n";
	print STDERR $cmd;`$cmd`;
	print STDERR "adjusted abundance file built\n";
}

####### compress the results to save space ### 
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

