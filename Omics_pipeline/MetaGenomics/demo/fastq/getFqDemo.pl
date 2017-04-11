#!/usr/bing/perl -w
use strict;

# got demo data for test
if(@ARGV<2){
	print "usage:\nperl $0 <gene_id_list> <sample_path_list(ID|soapPath|fqPath)>\n";
	exit;
}
my($gene_id_list,$sample_path_list) = @ARGV;
my(%INFO,@genes,@samples);
open ID,"< $gene_id_list";
open IF,"< $sample_path_list";
while(<IF>){
	chomp;
	next if $_ =~ /^#/;
	my @a = split;
	push @samples,$a[0];
	$INFO{$a[0]}{'soap'} = $a[1];
	$INFO{$a[0]}{'fastq'} = $a[2];
}
close IF;
while(<ID>){
	chomp;
	push @genes,$_;
}
close ID;
for(0..$#samples){
	my $s=$samples[$_];
	next if $s eq "";
	open SF,"gzip -dc $INFO{$s}{'soap'}|" or die $!;
	open RF,"gzip -dc $INFO{$s}{'fastq'}|" or die $!;
	open OT,"|gzip > $s.demo.fq.gz" or die $!;

#init
	my (@rawInfo, @soapInfo);
	my ($hitNum,$missNum,$c,$soapID,$hitGene) = (0,0,0,"","");
	while($c==0){
		@soapInfo   = split(/\t|\s/,<SF>);
		($soapID,$hitGene) = (&trimID($soapInfo[0]),$soapInfo[7]);
		$c = &ifSelect($hitGene);
	}

	while(<RF>){
		@rawInfo = ();
		push @rawInfo,$_;
		for(1..3){ $_ = <RF>; push @rawInfo,$_ }
		my $rawID = &trimID($rawInfo[0]);
		if($c==1 && $rawID eq $soapID && $hitNum < 10000){
			print OT @rawInfo;
			$hitNum ++;
			my $soapLine = "";
			$c = 0;
			while($c==0){
				for(1..$soapInfo[3]){$soapLine=<SF>};
				@soapInfo   = split(/\t|\s/,$soapLine);
				($soapID,$hitGene) = (&trimID($soapInfo[0]),$soapInfo[7]);
				$c = &ifSelect($hitGene);
			}
		}elsif($missNum < 10000){
			print OT @rawInfo;
			$missNum ++;
		}elsif($missNum >= 10000 && $hitNum >= 10000){
			last
		}
	}
	close RF;
	close SF;
	close OT;
}

##### SUB #####
sub trimID {
	$_ = shift;
	return("") if not defined $_;
	@_ = split /\t|\s/;
	$_[0] =~ s/^@//;
	$_[0] =~ s/\#(\w*(\/[12]|)|)$//;
	return($_[0]);
}

sub ifSelect {
	my $hitGene = shift;
	my $c = 0;
	for(0..$#genes){
		$c ++ ;
		last
	}
	return($c);
}

