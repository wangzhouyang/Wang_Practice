#!/usr/bin/perl
use strict;

my $usage = "USAGE:
	perl <map> <work dir> <steps> > outputfile
	map: seq_id\tsample_id.
	\@steps =(\"none\",\"filter\",\"rmhost\",\"soap\",\"abun\",\"soapALL\");\n";
die $usage if @ARGV <3;
my ($map,$wd,$step)=@ARGV;


$wd =~ s/\/$//;
my (%STAT,%MAP,%ABUN);
#directories location
my $d_clean = "$wd/clean";
my	$d_rmhost = "$wd/rmhost";
my	$d_soap = "$wd/soap";

my(@seqs,@inds)=();
open MAP,"<$map" or die $!;
while(<MAP>){
	chomp;
	my @a = split;
	$MAP{1}{$a[1]}=$a[0]; #individual -> sequence 
	$MAP{2}{$a[0]}=$a[1]; #sequence -> individual
	# In this script, @samples refers seq/indi depends on the objects under each steps.
	push @inds,$a[0];
	push @seqs,$a[1];
}

my @steps =("none","filter","rmhost","soap","abun","soapRepeat");
if($step =~ /1/){
	my @samples = @seqs;
	while ($#samples >-1){
		my(@heads,@vals)=();
		chomp(my $sam = shift @samples);
		open SC,"< $wd/clean/$sam.clean.stat_out" or die "can note open $wd/clean/$sam.clean.stat_out!".$!;
		chomp($_=<SC>);@heads = split /\t|\s/;
		chomp($_=<SC>);@vals = split /\t|\s/;
		@{$STAT{$sam}{1}{'H'}} = @heads;
		for (my $i=0;$i<=$#heads;$i++){
			$STAT{$sam}{1}{$heads[$i]} = $vals[$i];
		}
		close SC;
	}
}

if($step =~ /2/){
	my @samples = @seqs;
	while ($#samples >-1){
		my(@heads,@vals)=();
		chomp(my $sam = shift @samples);
		open SC,"< $wd/rmhost/$sam.rmhost.stat_out" or die $!; 
		chomp($_=<SC>);@heads = split /\t|\s/;
		chomp($_=<SC>);@vals = split /\t|\s/;
		@{$STAT{$sam}{2}{'H'}} = @heads;
		for (my $i=0;$i<=$#heads;$i++){
			$STAT{$sam}{2}{$heads[$i]} = $vals[$i];
		}
		close SC;
	}   
}

if($step =~ /3/){
	my @samples = @seqs;
	while ($#samples >-1){
		my(@heads,@vals)=();
		my $sam = shift @samples;
		my $path = "$wd/soap/$MAP{1}{$sam}.gene.build/$sam";
		my @logs = ();
		if (-e "$path.soap.SE.log"){push @logs,"SE"}
		if (-e "$path.soap.pair.log"){push @logs,"pair"}
		if (-e "$path.soap.single.log"){push @logs,"single"}
		while($#logs>-1){
			my $log = shift @logs;
            if($log eq "pair"){
                open SC,"tail -10 $path.soap.$log.log|" or die $!;
                @{$STAT{$sam}{"3.$log"}{'H'}} = ("$log.total","$log.paired","$log.singled");
                chomp($_=<SC>);$_ =~ /(\d+)/;
                $STAT{$sam}{"3.$log"}{"$log.total"} = $1;
                chomp($_=<SC>);$_ =~ /(\d+)/;
                $STAT{$sam}{"3.$log"}{"$log.paired"} = $1;
                chomp($_=<SC>);$_ =~ /(\d+)/;
                $STAT{$sam}{"3.$log"}{"$log.singled"} = $1;
                close SC;
            }else{
    			open SC,"tail -9 $path.soap.$log.log|" or die $!;
                @{$STAT{$sam}{"3.$log"}{'H'}} = ("$log.total","$log.aligned");
    			chomp($_=<SC>);$_ =~ /(\d+)/;
                $STAT{$sam}{"3.$log"}{"$log.total"} = $1;
                chomp($_=<SC>);$_ =~ /(\d+)/;
                $STAT{$sam}{"3.$log"}{"$log.aligned"} = $1;
			    close SC;
            }
		}
	}
}

if($step =~ /4/){
	my @samples = @inds;
	while ($#samples >-1){
		my $sam = shift @samples;
		chomp(my $size=`cat $wd/soap/$sam.abundance.size`);
		@{$STAT{$MAP{2}{$sam}}{4}{'H'}} = "size";
		$STAT{$MAP{2}{$sam}}{4}{'size'} = $size;
	}   
}

if($step =~ /A/){
	my @samples = @seqs;
	while ($#samples >-1){
		my(@heads,@vals)=();
		my $sam = shift @samples;
		my $path = "$wd/soapA/$MAP{1}{$sam}.gene.build/$sam";
		my @logs = ();
		if (-e "$path.soap.SE.log"){push @logs,"SE"}
		if (-e "$path.soap.pair.log"){push @logs,"pair"}
		if (-e "$path.soap.single.log"){push @logs,"single"}
		while($#logs>-1){
			my $log = shift @logs;
			open SC,"tail -9 $path.soap.$log.log|" or die $!;               
			chomp($_=<SC>);@heads = split /\s+|\t/;
			@{$STAT{$sam}{"5.$log"}{'H'}} = ("$log.reads","$log.aligned");
			$STAT{$sam}{"5.$log"}{"$log.reads"} = $heads[2];
			chomp($_=<SC>);@vals = split /\s+|\t/;
			$STAT{$sam}{"5.$log"}{"$log.aligned"} = $vals[1];
			close SC;
		}
	}
}  



my $time =1;
my $title = "id\tsample";
my $content = "";
foreach my $sam ( sort keys %STAT ){
	$content .="$sam\t$MAP{1}{$sam}";
	foreach my $part (sort keys %{$STAT{$sam}}){
		foreach my $head (@{$STAT{$sam}{$part}{'H'}}){
			#$label .= "\t$part" if $time;
			$title .= "\t$steps[$part]\_$head" if $time;
			$content .= (defined $STAT{$sam}{$part}{$head})?"\t$STAT{$sam}{$part}{$head}":"\tNA";
		}
	}
	$title .= "\n" if $time;
	$time = 0;
	$content .= "\n";
}
print "$title"."$content";









