#! /usr/bin/perl

=head1 Program: get_profile.pl

=head1 Version: V1.0

=head1 Updated: May 27 2013,  15:50:08

=head1 Description: This program use to get other profile from gene profile and annotation

=head1 
        	
	Usage: perl get_profile.pl [options]

        Options: 
        -g  <str>   parallelism.[first col: gene id; second col:annotation name]
	-f  <str>   gene profile file
	-o  <str>   output pfx.[report]

	Contact: liyin@genomics.org.cn
   
=head1 
     	
=cut

use strict;
use warnings;

use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
use Data::Dumper;
use Pod::Text;
#use threads;
#use threads::shared;

##initialize some parameters fot GetOptions
our ($parallelism, $gene_profile, $out, $anno_profile);
our (%parallelism_anno_hash, @line_ary, @sample_ary);
our (%anno_profile_hash, @anno_profile_ary, $i, $each);

GetOptions( 
        "g=s"=>\$parallelism,
	"f=s"=>\$gene_profile,
	"o=s"=>\$out,

);
##get the introduction information
die `pod2text $0` if ( !$parallelism || !$gene_profile);

#$out ||= ".";
#$out =~ s/\/$//;

my $pwd = $ENV{'PWD'};

#$out = "$pwd/$out" if($out !~ /^\//);

$parallelism = "$pwd/$parallelism" if($parallelism !~ /^\//);

$gene_profile = "$pwd/$gene_profile" if($gene_profile !~ /^\//);


#unless(-e $out){
#	`mkdir -p $out`;
#}

#profile result file
$anno_profile = $out;

#遍历获取每个基因ID对应的注释
open GS, $parallelism or die "can't read $parallelism:$!\n";
while(<GS>){
	chomp;
	@line_ary = split /\t/;
	$parallelism_anno_hash{$line_ary[0]} = $line_ary[1];  #记录gene id对应的注释

}
close GS;

#遍历gene profile
open GP, $gene_profile or die "can't read $gene_profile:$!\n";
$_=<GP>;
chomp;
@sample_ary = split(/\t/, $_);

while(<GP>){
	chomp;
	@line_ary = split /\t/;
	for($i = 1; $i<=$#line_ary; ++$i){  
		if($parallelism_anno_hash{$line_ary[0]}){
			$anno_profile_hash{$parallelism_anno_hash{$line_ary[0]}}{$sample_ary[$i]} += $line_ary[$i];  #统计每个注释每个样品的丰度
		}else{
			#$anno_profile_hash{"other"}[$sample_ary[$i]] += $line_ary[$i];          
		}
	}
}
close GP;


##输出profile
open annoP, ">$anno_profile" or die "can't write $anno_profile:$!\n";
for($i=0; $i<=$#sample_ary; ++$i){   #输出sample id到profile第一列
	print annoP $sample_ary[$i];
	if($i != $#sample_ary){
		print annoP "\t";
	}
}
print annoP "\n";

shift @sample_ary;
foreach my $anno (sort keys %anno_profile_hash){ #遍历输出anno profile
	my $anno_nospace = $anno;
	$anno_nospace =~ s/\s+/\_/g;
	print annoP $anno_nospace."\t";
	for( $i=0; $i<=$#sample_ary; ++$i){
		print annoP $anno_profile_hash{"$anno"}{$sample_ary[$i]};
		if($i != $#sample_ary){
			print annoP "\t";
		}
	}
	print annoP "\n";
}
close annoP;





