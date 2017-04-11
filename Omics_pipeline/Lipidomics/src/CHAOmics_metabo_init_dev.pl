#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use lib $Bin;
use Cwd 'abs_path';

sub usage {
	print <<USAGE
Version:	DEMO
Update :	20161008 	(since 20151120)
----------------------------------------
Description:
	This pipeline is to ease my pressure for Multiple omics analysis. For this version, I'm focus on the main analysis for publish.
Usage:
	perl $0 -cfg <configure> -dir <work dir> -step [options]
Options
	--help		ask for help
	--cfg	[essential] the configure file contains various fussy settings. Be patient.
	--dir	[essential] The shell and output date will be saved under this directory (default "project")
	--func	First 2 functions should be ran step by step, others are based on the scaled profile data.
			Do not seperate them by space or tabs.
			pretreat		formating dataset
			scale		use range-scale method to adjust data
			rank-sum			
			permanova	
			randomforest
e.g.
	perl $0 -cfg configureA -dir ./output -step polish

USAGE
}

# Globle variables
my ($help, $cfg, $dir, $step, $func);
GetOptions(
	"help"   => \$help,
	"cfg=s"  => \$cfg,
	"dir=s"  => \$dir,
	"step=s" => \$step,
	"func=s" => \$func,
);
#&usage && exit 0 if defined $help;
if ($help || not defined $cfg){
	&usage && exit 0;
}
# Default settings
$dir  ||= "project";
$dir = abs_path($dir);
# load tools & software configure
my (%library,%path);

# load config associated with the input data
my (%cfg,%dcfg);
my ($s_name,$ion_name,$position)=("","",0);
open IN,$cfg or die "fail open: $cfg\n";
while(<IN>){
	chomp;
	if (/###/){$position = 0 && next};
	if (/(\S+)\s*=\s*(\S+)/){
		my ($cfg_name,$cfg_value) = ($1,$2);
		if ($cfg_name eq "sample_pfx"){
			$s_name = $cfg_value;
			$ion_name = "";
			$position = 0;
		}elsif($cfg_name eq "ion"){
			$ion_name = $cfg_value;
			$position = 1;
		}elsif($position == 1){
			$cfg{$s_name}{$ion_name}{$cfg_name} = ($cfg_name =~ /prf|phe|lst|csv/)?abs_path($cfg_value):$cfg_value;
		}else {
			$dcfg{$cfg_name} = ($cfg_name =~ /phenotype|arg\.?\.phe/)?abs_path($cfg_value):$cfg_value;
		}
	}
}
close IN;

my $idmap_pl = "$Bin/idmapping.pl";
my $replace_pl = "$Bin/replace.pl";
my $trim_pl = "$Bin/trim.pl";
my $sort_pl = "$Bin/sortTable.pl";
my $Rscript = ($dcfg{'Rscript'})?$dcfg{'Rscript'}:"Rscript";

# setup directory
my $shdir=$dir.'/Shell';
system "mkdir -p $shdir" unless(-d $shdir);
my $batdir=$shdir.'/Batch';
system "mkdir -p $batdir" unless(-d $batdir);
my $samdir=$shdir.'/Sample';
system "mkdir -p $samdir " unless(-d $samdir);
#my $resdir=$dir.'/Result';
#system "mkdir -p $resdir" unless(-d $resdir);
## setup sub dir for each samples
foreach my $sample (sort keys %cfg) {
	my $specific_Samdir = $samdir."/".$sample;
	system "mkdir -p $specific_Samdir" unless(-d $specific_Samdir);
	foreach my $ion (sort keys %{$cfg{$sample}}){
		$cfg{$sample}{$ion}{'ion_dir'} = $specific_Samdir."/".$ion;
		system "mkdir -p $cfg{$sample}{$ion}{'ion_dir'}" unless(-d $cfg{$sample}{$ion}{'ion_dir'});
	}
}

######################################
#            FUNCTION PART           #
#          #################         #
(my $pfx_cfg = $cfg) =~ s/\.cfg//;
open FIN,">$dir/RUN.$pfx_cfg.sh";

# pretreatment of data source
#######################################
if ($func =~ /pretreat/) {
	my $batch_sh = $batdir."/00.pretreatment.$pfx_cfg.sh";
	my $result_dir = $dir."/00.data";
	system "mkdir -p $result_dir " unless(-d $result_dir);

	open BATCH,">$batch_sh";
	foreach my $sample (sort keys %cfg) {
		my $sample_sh = $samdir."/".$sample."/$sample.sh";
		my $unify_head_cmd = "\n";
		my $unify_tail_cmd = "\n";
		open UNIFY,">$sample_sh" or die "can't access to $sample_sh";
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $metabo_csv = $cfg{$sample}{$ion}{'metabo_csv'};
			my $new_ID_lst = $cfg{$sample}{$ion}{'new_ID_lst'};
			my $raw_ID_lst = $cfg{$sample}{$ion}{'raw_ID_lst'};
			my $ID_pfx     = $cfg{$sample}{$ion}{'ID_pfx'};
			my $ion_dir    = $cfg{$sample}{$ion}{'ion_dir'};
			my $res_pfx    = $result_dir."/".$sample."_".$ion;

			$cfg{$sample}{$ion}{'out.0.prf'} = "$res_pfx.pretreatment.res.comm.xls";

			my $pretreatment_sh = $ion_dir."/00.pretreatment.$sample.$ion.sh";
			open PT,">$pretreatment_sh" or die "can't access to $pretreatment_sh.$!\n";
			print PT "sh $Bin/S01_1.pretreatment.sh $metabo_csv $new_ID_lst $raw_ID_lst $res_pfx $ID_pfx $idmap_pl $replace_pl $trim_pl $sort_pl > $res_pfx.pretreatment.res.xls\n";
			close PT;
			$unify_head_cmd .= "head -1 $res_pfx.pretreatment.res.xls|xargs -n1 > $res_pfx.newID.lst\n$ion=$res_pfx.newID.lst\n";
			$unify_tail_cmd .= "perl $trim_pl $result_dir/$sample.comm.lst $res_pfx.pretreatment.res.xls in $cfg{$sample}{$ion}{'out.0.prf'} h &\n";
			print UNIFY "sh $pretreatment_sh &\n";
		}

		$dcfg{'out.0.phe'} = "$result_dir/$sample.comm.phenotype.tab";

		print UNIFY "wait\n$unify_head_cmd\ncomm \$pos \$neg -1 -2 > $result_dir/$sample.comm.lst\n$unify_tail_cmd\n";
		print UNIFY "perl $trim_pl $result_dir/$sample.comm.lst $dcfg{'phenotype'} in $dcfg{'out.0.phe'} r &\n";

		print UNIFY "wait\n";
		close UNIFY;
		print BATCH "sh $sample_sh\n";
	}
	close BATCH;
	print FIN "sh $batch_sh\n";
}

# scaling profile
###########################
if ($func =~ /scale/) {
	my $scaling_R = $Bin."/range_scaling.R";
	my $batch_sh = $batdir."/01.scaling.$pfx_cfg.sh";
	my $result_dir = $dir."/01.scaled_profile";
    system "mkdir -p $result_dir " unless(-d $result_dir);

	open BATCH,">$batch_sh";
	foreach my $sample (sort keys %cfg) {
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $ion_dir = $cfg{$sample}{$ion}{'ion_dir'};
			my $res_pfx = $result_dir."/".$sample."_".$ion;

			$cfg{$sample}{$ion}{'arg.1.prf'} ||= ($cfg{$sample}{$ion}{'out.0.prf'})?$cfg{$sample}{$ion}{'out.0.prf'}:"$dir/00.data/$sample.comm.phenotype.tab";
			$cfg{$sample}{$ion}{'out.1.prf'} = "$res_pfx.pretreatment.res.comm.range_scaling.xls";

			my $step_sh = $ion_dir."/01.scaling.$sample.$ion.sh";
			open STEP,">$step_sh" or die "can't access to $step_sh";
			print STEP "$Rscript $scaling_R $cfg{$sample}{$ion}{'arg.1.prf'} $cfg{$sample}{$ion}{'out.1.prf'}\n";
			print BATCH "sh $step_sh &\n";
			close STEP;
		}
	}
	print BATCH "wait\n";
	close BATCH;
	print FIN "sh $batch_sh\n";
}

# Rank-Sum test
##################
if ($func =~ /rank-sum/) {
	my $wilcox_R= "$Bin/wilcox.R";
	my $kw_R = "$Bin/Kruskal-Wallis.R";
    my $kw_direction_P = "$Bin/Kruskal-Wallis.direction.pl";
	my $batch_sh = $batdir."/02.wilcox.$pfx_cfg.sh";
	my $result_dir = $dir."/02.wilcox";
    system "mkdir -p $result_dir " unless(-d $result_dir);

	open BATCH,">$batch_sh";
	foreach my $sample (sort keys %cfg) {
		$dcfg{'arg.2.phe'} ||= ($dcfg{'out.0.phe'})?$dcfg{'out.0.phe'}:"$dir/00.data/$sample.comm.phenotype.tab";
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $ion_dir = $cfg{$sample}{$ion}{'ion_dir'};

			$dcfg{'arg.2.prf'} ||= ($cfg{$sample}{$ion}{'out.1.prf'})?$cfg{$sample}{$ion}{'out.1.prf'}:"$dir/01.scaled_profile/$sample\_$ion\.pretreatment.res.comm.range_scaling.xls"; #$dir."/01.scaled_profile/".$sample."_".$ion.".pretreatment.res.comm.range_scaling.xls";

			my @todo_phes  = split( /;/,$dcfg{'wilcox_todo_phe'});
			while (@todo_phes >0){
				my @cur_phe = split (/\|/,shift @todo_phes);
				
				my $todo_todo_phe = shift @cur_phe;
				my $todo_adjust = shift @cur_phe;
				
				my @todo_pairs = split( /;/,$dcfg{'wilcox_todo_pair'});
				while(@todo_pairs>0){
					my @todos = split(/,/, shift @todo_pairs);
#					(shift @todo_pairs) =~ /(\S+),(\S+)/;
					my $res_pfx = $result_dir."/".$sample."_".$ion;
					if (@todos == 2){
						my($todo_case, $todo_control) = ($todos[0], $todos[1]);
					
						my $step_sh = $ion_dir."/02.wilcox.$todo_todo_phe.$todo_case-vs-$todo_control.$sample.$ion.sh";
						open STEP,">$step_sh" or die "can't access to $step_sh";
						print STEP "$Rscript $wilcox_R $dcfg{'arg.2.prf'} $dcfg{'arg.2.phe'} $todo_case $todo_control $todo_todo_phe $res_pfx.$todo_todo_phe.$todo_case-vs-$todo_control.wilcox.xls\n";
						print BATCH "sh $step_sh &\n";
						close STEP;
					}elsif(@todos == 3){
						my($todo_ctr, $todo_case1, $todo_case2) = ($todos[0], $todos[1], $todos[2]);
						my $step_sh = $ion_dir."/02.Kw.$todo_todo_phe.$todo_ctr-vs-$todo_case1-vs-$todo_case2.adj_$todo_adjust.$sample.$ion.sh";
						open STEP,">$step_sh" or die "can't access to $step_sh\n$!\n";
						print STEP "$Rscript $kw_R $dcfg{'arg.2.prf'} $dcfg{'arg.2.phe'} $todo_ctr $todo_case1 $todo_case2 $todo_todo_phe $todo_adjust $res_pfx.$todo_todo_phe.$todo_ctr-vs-$todo_case1-vs-$todo_case2.adj_$todo_adjust\_kw.xls 1>$res_pfx.$todo_todo_phe.$todo_ctr-vs-$todo_case1-vs-$todo_case2.adj_$todo_adjust\_kw.log \n";
						print STEP "perl $kw_direction_P $res_pfx.$todo_todo_phe.$todo_ctr-vs-$todo_case1-vs-$todo_case2.adj_$todo_adjust\_kw.xls $res_pfx.$todo_todo_phe.$todo_ctr-vs-$todo_case1-vs-$todo_case2.adj_$todo_adjust\_kw.direction.xls\n";
						print BATCH "sh $step_sh &\n";
						close STEP;
					}
				}
			}
		}
	}
	close BATCH;
	print FIN "sh $batch_sh\n";
}


# permanova
# ###################
if ($func =~ /permanova/) {
	my $jsd_exc  = "$Bin/JSD_distance";
	my $bray_exc = "$Bin/euclidean_bray.distance";
	my $permanova_R = "$Bin/permanova.R";
	my $batch_sh = $batdir."/03.permanova.$pfx_cfg.sh";
	my $result_dir = $dir."/03.permanova";
	system "mkdir -p $result_dir " unless(-d $result_dir);
	my $dis_dir = $result_dir."/distance";
	system "mkdir -p $dis_dir  " unless(-d $dis_dir);

	open BATCH,">$batch_sh";
	foreach my $sample (sort keys %cfg) {
		$dcfg{'arg.3.phe'} ||= ($dcfg{'out.0.phe'})?$dcfg{'out.0.phe'}:"$dir/00.data/$sample.comm.phenotype.tab"; 
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $ion_dir = $cfg{$sample}{$ion}{'ion_dir'};
			
			my $dis_pfx = $dis_dir."/".$sample."_".$ion;
			my $res_pfx = $result_dir."/".$sample."_".$ion;
			my $step_sh = $ion_dir."/03.distance.$sample.$ion.sh";

			$cfg{$sample}{$ion}{'arg.3.prf'} ||= ($cfg{$sample}{$ion}{'out.1.prf'})?$cfg{$sample}{$ion}{'out.1.prf'}:"$dir/01.scaled_profile/$sample\_$ion.pretreatment.res.comm.range_scaling.xls";
			my $arg_distance_jsd_out = "$dis_pfx.JSD_distance";
			my $arg_distance_bray_out = "$dis_pfx.bray_distance";

			open STEP,">$step_sh" or die "can't access to $step_sh";
			print STEP "### distance building\n$jsd_exc $cfg{$sample}{$ion}{'arg.3.prf'} $arg_distance_jsd_out &\n";
			print STEP "$bray_exc $cfg{$sample}{$ion}{'arg.3.prf'} $dis_pfx &\n";
			close STEP;

			my $step2_sh = $ion_dir."/03.permanova.$sample.$ion.sh"; 
			open STEP2,">$step2_sh" or die "can't access to $step_sh";
			print STEP2 "wait\n### permanova\n";
			print STEP2 "$Rscript $permanova_R $arg_distance_jsd_out $dcfg{'arg.3.phe'} 5 $res_pfx.JSD &\n";
			print STEP2 "$Rscript $permanova_R $arg_distance_bray_out $dcfg{'arg.3.phe'} 5 $res_pfx.bary &\n";
			close STEP2;
			print BATCH "sh $step_sh\nsh $step2_sh &\n";
		}
	}
	close BATCH;
	print FIN "sh $batch_sh\n";
}

# randomForest
########################
if ($func =~ /randomforest/) {
	my $RF_ini_pl = "$Bin/randomForest_initial.pl";
	my $batch_sh = $batdir."/04.randomForest.$pfx_cfg.sh";
	my $result_dir = $dir."/04.randomForest";
	system "mkdir -p $result_dir " unless(-d $result_dir);

	open BATCH,">$batch_sh";
	foreach my $sample (sort keys %cfg) {
		$dcfg{'arg.4.phe'} ||= ($dcfg{'out.0.phe'})?$dcfg{'out.0.phe'}:"$dir/00.data/$sample.comm.phenotype.tab";
		my $arg_todo_phe = $dcfg{'RF_todo_phe'};
		$dcfg{'RF_todo_id'} ||= "F";
		$dcfg{'RF_todo_marker'} ||= "F";
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $ion_dir = $cfg{$sample}{$ion}{'ion_dir'};
			$cfg{$sample}{$ion}{'arg.4.prf'} ||= ($cfg{$sample}{$ion}{'out.1.prf'})?$cfg{$sample}{$ion}{'out.1.prf'}:"$dir/01.scaled_profile/$sample\_$ion.pretreatment.res.comm.range_scaling.xls";
			my $arg_ctr = $dcfg{'RF_todo_ctr'};
			my $arg_case1 = $dcfg{'RF_todo_case1'};
			my $arg_case2 = $dcfg{'RF_todo_case2'};
			my $arg_repeat = $dcfg{'RF_todo_repeat'};
			my $res_pfx = $result_dir."/".$sample."_".$ion;

			my $step_sh = $ion_dir."/04.randomForest.$sample.$ion.sh";
			#build scripts
			system("perl $RF_ini_pl -s $arg_repeat -phe $dcfg{'arg.4.phe'} -do_w T -todo $arg_todo_phe -pro $cfg{$sample}{$ion}{'arg.4.prf'} -todo_id $dcfg{'RF_todo_id'} -todo_marker $dcfg{'RF_todo_marker'} -ctr $arg_ctr -c1 $arg_case1 -c2 $arg_case2 -f 10 -ct 5 -cs 0.9 -pf $res_pfx");

			open STEP,">$step_sh" or die "can't access to $step_sh";

			print STEP "sh $res_pfx.batch/step1_setSamples.sh\n";
			print STEP "sh $res_pfx.batch/do_wilcox.sh\n";
			print STEP "sh $res_pfx.batch/step2_randomForest.sh\n";
			close STEP;
			print BATCH "sh $step_sh &\n";
		}
	}
	close BATCH;
	print FIN "sh $batch_sh\n";
}

# correlation
########################
if ($func =~ /correlation/){
	my $COR_R = "$Bin/docorelation.R";
	my $batch_sh = $batdir."/05.correlation.$pfx_cfg.sh";
	my $result_dir = $dir."/05.correlation";
	system "mkdir -p $result_dir " unless(-d $result_dir);

	open BATCH,">$batch_sh";
	foreach my $sample (sort keys %cfg){
		$dcfg{'arg.5.phe'} ||=($dcfg{'out.0.phe'})?$dcfg{'out.0.phe'}:"$dir/00.data/$sample.comm.phenotype.tab";
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $ion_dir = $cfg{$sample}{$ion}{'ion_dir'};
			$cfg{$sample}{$ion}{'arg.5.prf'} ||= ($cfg{$sample}{$ion}{'out.1.prf'})?$cfg{$sample}{$ion}{'out.1.prf'}:"$dir/01.scaled_profile/$sample\_$ion.pretreatment.res.comm.range_scaling.xls";
			my $res_pfx = $result_dir."/".$sample."_".$ion;

			my $step_sh = $ion_dir."/05.correlation.$sample.$ion.sh";
			#build scripts
			open STEP,">$step_sh" or die "can't access to $step_sh\n";
			print STEP "$Rscript $COR_R $dcfg{'arg.5.phe'} $cfg{$sample}{$ion}{'arg.5.prf'} $res_pfx.spearman.xls\n";
			close STEP;
			print BATCH "sh $step_sh\n";
		}
	}
	close BATCH;
	print FIN "sh $batch_sh\n";
}

#glm
if ($func =~ /glm/){
	my $GLM_R = "$Bin/GLM.logistic.R";
	my $batch_sh = $batdir."/06.GLM.$pfx_cfg.sh";
	my $result_dir = $dir."/06.GLM";
	system "mkdir -p $result_dir " unless(-d $result_dir);

	open BATCH,">$batch_sh";
	foreach my $sample (sort keys %cfg){
		$dcfg{'arg.6.phe'} ||=($dcfg{'out.0.phe'})?$dcfg{'out.0.phe'}:"$dir/00.data/$sample.comm.phenotype.tab";
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $ion_dir = $cfg{$sample}{$ion}{'ion_dir'};
			$cfg{$sample}{$ion}{'arg.6.prf'} ||= ($cfg{$sample}{$ion}{'out.1.prf'})?$cfg{$sample}{$ion}{'out.1.prf'}:"$dir/01.scaled_profile/$sample\_$ion.pretreatment.res.comm.range_scaling.xls";
			my $res_pfx = $result_dir."/".$sample."_".$ion;

			my $step_sh = $ion_dir."/06.GLM.$sample.$ion.sh";
			#build scripts
			open STEP,">$step_sh" or die "can't access to $step_sh\n";
			print STEP "$Rscript $GLM_R $cfg{$sample}{$ion}{'arg.6.prf'} $dcfg{'arg.6.phe'} $res_pfx\n";
			close STEP;
			print BATCH "sh $step_sh\n";
		}
	}
	close BATCH;
	print FIN "sh $batch_sh\n";
}

##########
close FIN;





######################################
#          SUB FUNCTION PART         #
#          #################         #

#usage: $parse_cfg($cfg_file,$cfg_hash,'circulation key variable'(if need));
sub read_cfg {
	my ($cfg_file, $cfg_hash, $key) = @_;
	open IN,"<$cfg_file" or die "fail open: $cfg_file";
	while(<IN>){
		chomp;
		next if (/^#/);
		if (/(\S+)\s*=\s*(\S+)/){
			my ($cfg_name,$cfg_value) = ($1,$2);
			$cfg_hash -> {$cfg_name} = $cfg_value;
		}
	}
	close IN;
}

#####################################################################
# update info
# 
# 20151202      complete steps: inputs pretreatments and formating
#                               scaling
#                               wilcox
#                               permanova
#                               randomeforest
#20151120      start workframe building
