#!/usr/bin/perl

=head1
	
	TBseq - a computational pipeline for detecting variants in NGS-data

	Copyright (C) 2016 Thomas A. Kohl, Maria R. De Filippo, Robin Koch, Viola Schleusener, Christian Utpatel, Daniela M. Cirillo, Stefan Niemann

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut


use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use Cwd;
use Getopt::Long;
use TBreads;
use TBbwa;
use TBstats;
use TBmerge;
use TBrefine;
use TBpile;
use TBtable;
use TBvariants;
use TBjoin;
use TBamend;
use TBstrains;
use TBgroups;
use TBtools;

# Specials for MTBC.
my  	$res             	=	"$RealBin/var/res/Base_Calibration_List.vcf"; # Known SNP positions for base call recalibration.
my	$resi_list_master	=	"$RealBin/var/res/MTB_Resistance_Mediating.txt"; # MTBC Resistance mediating variants and MTBC phylogentic SNPs.
my	$int_regions		=	"$RealBin/var/res/MTB_Extended_Resistance_Mediating.txt"; # MTBC extended intergenic resistance mediating positions.
my	$categories		=	"$RealBin/var/cat/MTB_Gene_Categories.txt"; # MTBC essential and non-essential genes.

# Get current working directory and time.
my  	$W_dir			=	getcwd();
my	$date_string		=	timer();
$date_string                   	=~      s/\[//;
$date_string                    =~      s/\]//;
$date_string                    =~      s/\-/_/g;
$date_string                    =~      s/\:/_/g;
$date_string                    =~      s/\s.*//;
$date_string                    =~      s/_/-/g;

# Define directories. Will be exported to the certain modules for easier updating purposes.
my	$BAM_OUT         	=	"$W_dir/Bam";
my 	$MBAM_OUT        	=	"$W_dir/Merged_Bam";
my 	$GATK_OUT        	=	"$W_dir/GATK_Bam";
my 	$MPILE_OUT       	=	"$W_dir/Mpileup";
my 	$POS_OUT         	=	"$W_dir/Position_Tables";
my 	$CALL_OUT        	=	"$W_dir/Called";
my	$JOIN_OUT		=	"$W_dir/Joint";
my	$AMEND_OUT		=	"$W_dir/Amend";
my	$STRAIN_OUT		=	"$W_dir/Classification";
my	$GROUPS_OUT		=	"$W_dir/Groups";

# Define PATHs. Will be exported to the certain modules for easier updating purposes.
my  	$VAR_dir		=	"$RealBin/var/ref";
my  	$BWA_dir         	= 	"$RealBin/opt/bwa_0.7.15";
my	$SAMTOOLS_dir    	=	"$RealBin/opt/samtools_1.3.1";
my	$PICARD_dir		=	"$RealBin/opt/picard_2.7.1";
my	$SAMBAMBA_dir    	=	"$RealBin/opt/sambamba_0.6.5";
my	$IGV_dir         	=	"$RealBin/opt/IGVTools_2.3.88";
my 	$GATK_dir        	=	"$RealBin/opt/GenomeAnalysisTK_3.6";

# Initialize commandline parameter.
my  	$step			=	"";
my  	$continue		=	"";
my      $samples                =       "";
my      $group_name             =       "";
my  	$ref			=	"";
my  	$machine		=	"";
my  	$run_number		=	"";
my  	$seq_length		=	"";
my  	$output_mode		=	"";
my	$mibqual		=	"";
my	$micovf			=	"";
my	$micovr			=	"";
my 	$miphred20		=	"";
my	$mifreq			=	"";
my	$unambigous		=	"";
my	$window			=	"";
my	$distance		=	"";
my  	$threads		=	"";
my	$lowmem			=	"";
my  	$help			=	"";

# Get the commandline parameter.
GetOptions	(	'step:s'        =>	\$step,
			'samples:s'	=>	\$samples,
			'groupname:s'	=>	\$group_name,
			'continue'      =>	\$continue,
			'ref:s'         =>	\$ref,
			'machine:s'     =>	\$machine,
			'run:s'         =>	\$run_number,
			'readlen:i'     =>	\$seq_length,
			'outmode:i'     =>	\$output_mode,
			'minbqual:i'	=>	\$mibqual,
			'mincovf:i'	=>	\$micovf,
			'mincovr:i'	=>	\$micovr,
			'minphred20:i'	=>	\$miphred20,
			'minfreq:i'	=>	\$mifreq,
			'unambig:i'	=>	\$unambigous,
			'window:i'	=>	\$window,
			'distance:i'	=>	\$distance,
			'threads:i'     =>	\$threads,
			'lowmem'	=>	\$lowmem,
            		'help'          =>	\$help
		);

# Print help message if specified or error if step is not defined or wrong.
if	($help eq '1'		)       	{	help();		exit 1;		}
if	($step eq ''		)       	{	nostep();	exit 1;		}
unless	(($step eq 'TBfull'	) 	||
	 ($step eq 'TBreads'	)	||
	 ($step eq 'TBbwa'	)	||
	 ($step eq 'TBstats' 	)	||
	 ($step eq 'TBmerge'	)	||
	 ($step eq 'TBrefine'	)	||
	 ($step eq 'TBpile'	)  	||
	 ($step eq 'TBtable'	)	||
	 ($step eq 'TBvariants'	)	||
	 ($step eq 'TBjoin'    	)	||
	 ($step eq 'TBamend'   	)	||
	 ($step eq 'TBstrains'  )	||
	 ($step eq 'TBgroups'	)
  	)					{	badstep($step);	exit 1;		}

# Test the commandline options and set default vaules if necessary.
if($continue		eq	''	)	{	$continue       =	0;					}
if($samples             eq      ''      )       {       $samples        =       "NONE";     				}
if($group_name          eq      ''      )       {       $group_name     =       "NONE";                           	}
if($ref			eq	''	)	{	$ref            =	"M._tuberculosis_H37Rv_2015-11-13";	}
if($machine		eq	''	)	{	$machine        =	"NGS-MACHINE";				}
if($run_number		eq	''	)	{	$run_number   	=	"nXXXX";				}
if($seq_length		eq 	''	)	{	$seq_length   	=	"151";					}
if($output_mode		eq 	''	)	{	$output_mode	=	0;					}
if($output_mode		>  	'3'	)	{	$output_mode    =	0;					}
if($mibqual		eq	''	)	{	$mibqual        =	13;					}
if($micovf		eq	''	)	{	$micovf         =	4;					}
if($micovr		eq	''	)	{	$micovr         =	4;                                      }
if($miphred20		eq	''	)	{	$miphred20	=	4;                                      }
if($mifreq		eq	''	)	{	$mifreq         =	75;					}
if($unambigous		eq	''	)	{	$unambigous	=	95;					}
if($window		eq	''	)	{	$window		=	12;					}
if($distance		eq	''	)	{	$distance	=	12;					}
if($lowmem		eq	''	)	{	$lowmem		=	0;					}
if($threads		eq	''	)	{	$threads	=	1;					}

# Set the name of the $ref fasta file and gene annotation.
my $refg     =      $ref;
$ref        .=      ".fasta";
$refg       .=      "_genes.txt";

# If $ref is not MTB than turn off resistance check.
unless($ref eq 'M._tuberculosis_H37Rv_2015-11-13.fasta') {
	$resi_list_master	=	"NONE";
	$int_regions		=	"NONE";
	$categories		=	"NONE";
}

# Inform the user what is used in general.
print  "\n<INFO>\t",timer(),"\tYou are $ENV{USER}.\n";
print  "<INFO>\t",timer(),"\tYour current working directory is: $W_dir\n";
print  "<INFO>\t",timer(),"\tYou requested $threads threads for the pipeline.\n";

print  "\n<INFO>\t",timer(),"\tYour parameter setting is:\n";
print  "<INFO>\t",timer(),"\t--step\t\t$step\n";
print  "<INFO>\t",timer(),"\t--continue\t$continue\n";
print  "<INFO>\t",timer(),"\t--samples\t$samples\n";
print  "<INFO>\t",timer(),"\t--groupname\t$group_name\n";
print  "<INFO>\t",timer(),"\t--ref\t\t$ref\n";
print  "<INFO>\t",timer(),"\t--machine\t$machine\n";
print  "<INFO>\t",timer(),"\t--run\t\t$run_number\n";
print  "<INFO>\t",timer(),"\t--readlen\t$seq_length\n";
print  "<INFO>\t",timer(),"\t--outmode\t$output_mode\n";
print  "<INFO>\t",timer(),"\t--minbqual\t$mibqual\n";
print  "<INFO>\t",timer(),"\t--mincovf\t$micovf\n";
print  "<INFO>\t",timer(),"\t--mincovr\t$micovr\n";
print  "<INFO>\t",timer(),"\t--minphred20\t$miphred20\n";
print  "<INFO>\t",timer(),"\t--minfreq\t$mifreq\n";
print  "<INFO>\t",timer(),"\t--unambig\t$unambigous\n";
print  "<INFO>\t",timer(),"\t--window\t$window\n";
print  "<INFO>\t",timer(),"\t--distance\t$distance\n";
print  "<INFO>\t",timer(),"\t--lowmem\t$lowmem\n";

print  "\n<INFO>\t",timer(),"\tFollowing programs will be used, if necessary:\n"; 
print  "<INFO>\t",timer(),"\t$BWA_dir\n";
print  "<INFO>\t",timer(),"\t$SAMTOOLS_dir\n";
print  "<INFO>\t",timer(),"\t$PICARD_dir\n";
print  "<INFO>\t",timer(),"\t$SAMBAMBA_dir\n";
print  "<INFO>\t",timer(),"\t$IGV_dir\n";
print  "<INFO>\t",timer(),"\t$GATK_dir\n";

print  "\n<INFO>\t",timer(),"\tFollowing directories will be used, if necessary:\n";
print  "<INFO>\t",timer(),"\t$BAM_OUT\n";
print  "<INFO>\t",timer(),"\t$MBAM_OUT\n";
print  "<INFO>\t",timer(),"\t$GATK_OUT\n";
print  "<INFO>\t",timer(),"\t$MPILE_OUT\n";
print  "<INFO>\t",timer(),"\t$POS_OUT\n";
print  "<INFO>\t",timer(),"\t$CALL_OUT\n";
print  "<INFO>\t",timer(),"\t$JOIN_OUT\n";
print  "<INFO>\t",timer(),"\t$AMEND_OUT\n";
print  "<INFO>\t",timer(),"\t$STRAIN_OUT\n";
print  "<INFO>\t",timer(),"\t$GROUPS_OUT\n";

system("mkdir -p $BAM_OUT");
system("mkdir -p $MBAM_OUT");
system("mkdir -p $GATK_OUT");
system("mkdir -p $MPILE_OUT");
system("mkdir -p $POS_OUT");
system("mkdir -p $CALL_OUT");
system("mkdir -p $JOIN_OUT");
system("mkdir -p $AMEND_OUT");
system("mkdir -p $STRAIN_OUT");
system("mkdir -p $GROUPS_OUT");

# Initialize content arrays.
my @fastq_files;
my @fastq_files_new;
my @bam_files;
my @gatk_files;
my @gatk_files_new;
my @mpile_files;
my @mpile_files_new;
my @pos_files;
my @pos_files_new;
my @var_files;
my @var_files_new;
my @join_files;
my @join_files_new;
my @amend_files;
my @amend_files_new;
my @group_files;
my @samples;
my %check_up;
if(-f "$W_dir/$samples") {
        open(IN,"$W_dir/$samples") || die "<ERROR>\t",timer(),"\tCan\'t open $samples: $!";
        @samples        =       <IN>;
        close(IN);
}
my $sample_number       =       scalar(@samples);

# Juming to certain pipeline steps.
if($step	eq	'TBfull'	)	{   goto    TBfull	}
if($step	eq	'TBreads'	)	{   goto    TBreads	}
if($step	eq	'TBbwa'		)	{   goto    TBbwa	}
if($step	eq	'TBmerge'	)	{   goto    TBmerge	}
if($step	eq	'TBrefine'	)	{   goto    TBrefine	}
if($step	eq	'TBpile'	)	{   goto    TBpile	}
if($step	eq	'TBtable'	)	{   goto    TBtable	}
if($step	eq	'TBvariants'	)	{   goto    TBvariants	}
if($step        eq      'TBstats'       )       {   goto    TBstats     }
if($step	eq	'TBjoin'	)	{   goto    TBjoin	}
if($step	eq	'TBamend'	)	{   goto    TBamend	}
if($step	eq	'TBstrains'	)	{   goto    TBstrains   }
if($step        eq      'TBgroups'      )       {   goto    TBgroups    }

# Pipeline execution.
TBfull:
print  "\n<INFO>\t",timer(),"\t### [TBfull] workflow selected! ###\n";



TBreads:
if($step eq 'TBreads') {
	print  "\n<INFO>\t",timer(),"\t### [TBreads] workflow selected ###\n";
}	
opendir(WORKDIR,"$W_dir") || die "<ERROR>\t",timer(),"\tCan\'t open directory $W_dir: $!";
@fastq_files		=	grep { $_ =~ /^\w.*\.fastq\.gz$/ && -f "$W_dir/$_" } readdir(WORKDIR);
if(scalar(@fastq_files) == 0) {
	print  "\n<ERROR>\t",timer(),"\tNo files to process! Check content of $W_dir!\n";
	exit 1;
}
print  "\n<INFO>\t",timer(),"\tProcessing read files:\n";
foreach my $fastq (sort { $a cmp $b } @fastq_files) {
	print  "<INFO>\t",timer(),"\t$fastq\n";
}
print  "\n<INFO>\t",timer(),"\tStart read file processing logic...\n";
tbreads($W_dir,$machine,$run_number,$seq_length,@fastq_files);
print  "<INFO>\t",timer(),"\tFinished read file processing logic!\n";
@fastq_files		=	();
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBbwa:
if($step eq 'TBbwa') {
	print  "\n<INFO>\t",timer(),"\t### [TBbwa] workflow selected ###\n";
}
opendir(WORKDIR,"$W_dir") || die "<ERROR>\t",timer(),"\tCan\'t open directory $W_dir: $!";
opendir(BAMDIR,"$BAM_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $BAM_OUT: $!";
@fastq_files		=	grep { $_ =~ /^\w.*\.fastq\.gz/ && -f "$W_dir/$_" } readdir(WORKDIR);
@bam_files              =       grep { $_ =~ /^\w.*\.bam$/ && -f "$BAM_OUT/$_" } readdir(BAMDIR);
if(scalar(@fastq_files) == 0) {
        print  "\n<ERROR>\t",timer(),"\tNo files to map! Check content of $W_dir!\n";
        exit 1;
}
%check_up		=	map { (my $id = $_) =~ s/\.bam$//; $id => $id; } @bam_files;
for(my $i = 0;$i < scalar(@fastq_files);$i++) {
	my $tmp		=	$fastq_files[$i];
	$tmp		=~	s/_R\d\.fastq.gz//;
	if(!exists $check_up{$tmp}) {
		push(@fastq_files_new,$fastq_files[$i]);
	}
}
if(scalar(@fastq_files_new) == 0) {
        print  "\n<ERROR>\t",timer(),"\tAll .fastq.gz files have been already mapped from $W_dir!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tMapping samples:\n";
foreach my $fastq (sort { $a cmp $b } @fastq_files_new) {
	print  "<INFO>\t",timer(),"\t$fastq\n";
}
print  "\n<INFO>\t",timer(),"\tStart BWA mapping logic...\n";
tbbwa($W_dir,$VAR_dir,$BWA_dir,$SAMTOOLS_dir,$BAM_OUT,$ref,$threads,@fastq_files_new);
print  "<INFO>\t",timer(),"\tFinished BWA mapping logic!\n";
@fastq_files		=	();
@fastq_files_new	=	();
@bam_files		=	();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBmerge:
if($step eq 'TBmerge') {
	print  "\n<INFO>\t",timer(),"\t### [TBmerge] workflow selected ###\n";
}
opendir(BAMDIR,"$BAM_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $BAM_OUT: $!";
@bam_files 		= 	grep { $_ =~ /^\w.*\.bam$/ && -f "$BAM_OUT/$_" } readdir(BAMDIR);
if(scalar(@bam_files) == 0) {
   	print  "\n<ERROR>\t",timer(),"\tNo files to merge! Check content of $BAM_OUT\n";
    	exit 1;
}
print  "\n<INFO>\t",timer(),"\tMerging mappings:\n";
foreach my $bam (sort { $a cmp $b } @bam_files) {
	print  "<INFO>\t",timer(),"\t$bam\n";
}
print  "\n<INFO>\t",timer(),"\tStart merging logic...\n";
tbmerge($SAMBAMBA_dir,$BAM_OUT,$MBAM_OUT,$GATK_OUT,$MPILE_OUT,$POS_OUT,$CALL_OUT,$threads,@bam_files);
print  "<INFO>\t",timer(),"\tFinished merging logic!\n";
@bam_files		=	();
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBrefine:
if($step eq 'TBrefine') {
	print  "\n<INFO>\t",timer(),"\t### [TBrefine] workflow selected ###\n";
}
opendir(BAMDIR,"$BAM_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $BAM_OUT: $!";
opendir(GATKDIR,"$GATK_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $GATK_OUT: $!";
@bam_files 		= 	grep { $_ =~ /^\w.*\.bam$/ && -f "$BAM_OUT/$_" } readdir(BAMDIR);
@gatk_files		=	grep { $_ =~ /^\w.*\.gatk\.bam$/ && -f "$GATK_OUT/$_" } readdir(GATKDIR);
if(scalar(@bam_files) == 0) {
    	print  "\n<ERROR>\t",timer(),"\tNo files to refine! Check content of $BAM_OUT!\n";
    	exit 1;
}
%check_up               =       map { (my $id = $_) =~ s/\.gatk\.bam$//; $id => $id; } @gatk_files;
for(my $i = 0;$i < scalar(@bam_files);$i++) {
        my $tmp         =       $bam_files[$i];
        $tmp            =~      s/\.bam//;
        if(!exists $check_up{$tmp}) {
                push(@gatk_files_new,$bam_files[$i]);
        }
}
if(scalar(@gatk_files_new) == 0) {
        print  "\n<ERROR>\t",timer(),"\tAll .bam files have been already refined from $BAM_OUT!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tRefining mappings:\n";
foreach my $bam (sort { $a cmp $b } @gatk_files_new) {
	print  "<INFO>\t",timer(),"\t$bam\n";
}
print  "\n<INFO>\t",timer(),"\tStart GATK refinement logic...\n";
tbrefine($W_dir,$VAR_dir,$PICARD_dir,$IGV_dir,$GATK_dir,$BAM_OUT,$GATK_OUT,$ref,$res,$threads,@gatk_files_new);
print  "<INFO>\t",timer(),"\tFinished GATK refinement logic!\n";
@bam_files		=	();
@gatk_files		=	();
@gatk_files_new		=	();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBpile:
if($step eq 'TBpile') {
	print  "\n<INFO>\t",timer(),"\t### [TBpile] workflow selected ###\n";
}
opendir(GATKDIR,"$GATK_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $GATK_OUT: $!";
opendir(MPILEDIR,"$MPILE_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $MPILE_OUT: $!";
@gatk_files 		= 	grep { $_ =~ /^\w.*\.gatk\.bam$/ && -f "$GATK_OUT/$_" } readdir(GATKDIR);
@mpile_files		=	grep { $_ =~ /^\w.*\.gatk\.mpileup$/ && -f "$MPILE_OUT/$_" } readdir(MPILEDIR);
if(scalar(@gatk_files) == 0) {
	print  "\n<ERROR>\t",timer(),"\tNo files to create mpileup files! Check content of $GATK_OUT!\n";
	exit 1;
}
%check_up               =       map { (my $id = $_) =~ s/\.mpileup$//; $id => $id; } @mpile_files;
for(my $i = 0;$i < scalar(@gatk_files);$i++) {
        my $tmp         =       $gatk_files[$i];
        $tmp            =~      s/\.bam//;
        if(!exists $check_up{$tmp}) {
                push(@gatk_files_new,$gatk_files[$i]);
        }
}
if(scalar(@gatk_files_new) == 0) {
        print  "\n<ERROR>\t",timer(),"\tAll .mpileup files have been already created from $GATK_OUT!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tCreating mpileups:\n";
foreach my $bam (sort { $a cmp $b } @gatk_files_new) {
	print  "<INFO>\t",timer(),"\t$bam\n";
}
print  "\n<INFO>\t",timer(),"\tStart creating mpileup files logic...\n";
tbpile($VAR_dir,$SAMTOOLS_dir,$GATK_dir,$GATK_OUT,$MPILE_OUT,$ref,$threads,@gatk_files_new);
print  "<INFO>\t",timer(),"\tFinished creating mpileup files logic!\n";
@gatk_files		=	();
@mpile_files		=	();
@gatk_files_new		=	();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBtable:
if($step eq 'TBtable') {
    	print  "\n<INFO>\t",timer(),"\t### [TBtable] workflow selected ###\n";
}
opendir(MPILEDIR,"$MPILE_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $MPILE_OUT: $!";
opendir(POSDIR,"$POS_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $POS_OUT: $!";
@mpile_files 		= 	grep { $_ =~ /^\w.*\.gatk\.mpileup$/ && -f "$MPILE_OUT/$_" } readdir(MPILEDIR);
@pos_files		=	grep { $_ =~ /^\w.*\.gatk_position_table\.tab$/ && -f "$POS_OUT/$_" } readdir(POSDIR);
if(scalar(@mpile_files) == 0) {
    	print  "\n<ERROR>\t",timer(),"\tNo files to create position tables! Check content of $MPILE_OUT!\n";
    	exit 1;
}
%check_up               =       map { (my $id = $_) =~ s/\.gatk_position_table\.tab$//; $id => $id; } @pos_files;
for(my $i = 0;$i < scalar(@mpile_files);$i++) {
        my $tmp         =       $mpile_files[$i];
        $tmp            =~      s/\.gatk\.mpileup//;
        if(!exists $check_up{$tmp}) {
                push(@mpile_files_new,$mpile_files[$i]);
        }
}
if(scalar(@mpile_files_new) == 0) {
        print  "\n<ERROR>\t",timer(),"\tAll position tables have been already created from $MPILE_OUT!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tCreating position tables:\n";
foreach my $mpile (sort { $a cmp $b } @mpile_files_new) {
	print  "<INFO>\t",timer(),"\t$mpile\n";
}
print  "\n<INFO>\t",timer(),"\tStart creating position tables logic...\n";
tbtable($VAR_dir,$MPILE_OUT,$POS_OUT,$ref,$mibqual,$lowmem,$threads,@mpile_files_new);
print  "<INFO>\t",timer(),"\tFinished creating position tables logic!\n";
@mpile_files		=	();
@pos_files		=	();
@mpile_files_new	=	();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBvariants:
if($step eq 'TBvariants') {
    	print  "\n<INFO>\t",timer(),"\t### [TBvariants] workflow selected ###\n";
}
opendir(POSDIR,"$POS_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $POS_OUT: $!";
opendir(CALLDIR,"$CALL_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $CALL_OUT: $!";
@pos_files 		= 	grep { $_ =~ /^\w.*\.gatk_position_table\.tab$/ && -f "$POS_OUT/$_" } readdir(POSDIR);
@var_files		=	grep { $_ =~ /^\w.*\.gatk_position_variants_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode$output_mode\.tab$/ && -f "$CALL_OUT/$_" } readdir(CALLDIR);
if(scalar(@pos_files) == 0) {
    	print  "\n<ERROR>\t",timer(),"\tNo files to call variants! Check content of $POS_OUT!\n";
    	exit 1;
}
%check_up               =       map { (my $id = $_) =~ s/\.gatk_position_variants_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode$output_mode\.tab$//; $id => $id; } @var_files;
for(my $i = 0;$i < scalar(@pos_files);$i++) {
        my $tmp         =       $pos_files[$i];
        $tmp            =~      s/\.gatk_position_table\.tab//;
        if(!exists $check_up{$tmp}) {
                push(@pos_files_new,$pos_files[$i]);
        }
}
if(scalar(@pos_files_new) == 0) {
        print  "\n<ERROR>\t",timer(),"\tAll variants files have been already created from $CALL_OUT!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tCalling variants:\n";
foreach my $pos_file (sort { $a cmp $b } @pos_files_new) {
	print  "<INFO>\t",timer(),"\t$pos_file\n";
}
print  "\n<INFO>\t",timer(),"\tStart variant calling logic...\n";
tbvariants($VAR_dir,$POS_OUT,$CALL_OUT,$ref,$refg,$micovf,$micovr,$miphred20,$mifreq,$resi_list_master,$int_regions,$output_mode,@pos_files_new);
print  "<INFO>\t",timer(),"\tFinished variant calling logic!\n";
@pos_files		=	();
@var_files		=	();
@pos_files_new		=	();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBstats:
if($step eq 'TBstats') {
        print  "\n<INFO>\t",timer(),"\t### [TBstats] workflow selected ###\n";
}
opendir(BAMDIR,"$BAM_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $BAM_OUT: $!";
@bam_files              =       grep { $_ =~ /^\w.*\.bam$/ && -f "$BAM_OUT/$_" } readdir(BAMDIR);

if(scalar(@bam_files) == 0) {
        print  "\n<ERROR>\t",timer(),"\tNo files to calculate statistics! Check content of $BAM_OUT!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tCalculating statistics for mappings:\n";
foreach my $bam (sort { $a cmp $b } @bam_files) {
        print  "<INFO>\t",timer(),"\t$bam\n";
}
print  "\n<INFO>\t",timer(),"\tStart statistics logic...\n";
tbstats($W_dir,$VAR_dir,$SAMTOOLS_dir,$BAM_OUT,$ref,$date_string,@bam_files);
print  "<INFO>\t",timer(),"\tFinished statistics logic!\n";
@bam_files              =       ();
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBjoin:
if($step eq 'TBjoin') {
    	print  "\n<INFO>\t",timer(),"\t### [TBjoin] workflow selected ###\n";
}
if($samples eq 'NONE') {
        print  "\n<INFO>\t",timer(),"\tSkipping Joint analysis. No samples defined for joint analysis.\n";
        print  "\n<INFO>\t",timer(),"\t### TBseq finished!!! ###\n";
        exit 1;
}
open(IN,"$W_dir/$samples") || die "<INFO>\t",timer(),"\tCan't find $samples file in $W_dir. Check $W_dir!\n";
while(<IN>) {
	my $line	=	$_;
	$line		=~	s/\015?\012?$//;
	my @fields	=	split(/\t/,$line);
	$check_up{$fields[0]."_".$fields[1]}	=   $fields[0]."_".$fields[1];
}
close(IN);
opendir(CALLDIR,"$CALL_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $CALL_OUT: $!";
opendir(JOINDIR,"$JOIN_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $JOIN_OUT: $!";
@var_files		=	grep { $_ =~ /^\w.*\.gatk_position_variants_cf$micovf\_cr$micovf\_fr$mifreq\_ph$miphred20\_outmode$output_mode\.tab$/ && -f "$CALL_OUT/$_" } readdir(CALLDIR);
@join_files		=	grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode3_samples$sample_number\.tab$/ && -f "$JOIN_OUT/$_" } readdir(JOINDIR);
if(scalar(@var_files) == 0) {
	print  "\n<ERROR>\t",timer(),"\tNo files to create joint variant tables! Check content of $POS_OUT and $CALL_OUT!\n";
	exit 1;
}
if(scalar(@join_files) != 0) {
	print  "\n<ERROR>\t",timer(),"\tIt seems that you already created a joint analysis for group $group_name and submitted parameter setting! Please check content of $JOIN_OUT!\n";
	exit 1;
}
for(my $i = 0;$i < scalar(@var_files);$i++) {
        my @tmp         =       split(/_/,$var_files[$i]);
        if(exists $check_up{$tmp[0]."_".$tmp[1]}) {
		print  "\n<INFO>\t",timer(),"\tFound $var_files[$i]! Using it for joint analysis...";
                push(@var_files_new,$var_files[$i]);
		delete($check_up{$tmp[0]."_".$tmp[1]});
        }
}
foreach my $samplelib(keys %check_up) {
	print  "\n<ERROR>\t",timer(),"\t Do you run already TBtable and TBvariants for $samplelib?\n";
	exit 1;
}
print  "\n\n<INFO>\t",timer(),"\tCreating a joint variant table from:\n";
foreach my $var_file (sort { $a cmp $b } @var_files_new) {
    	print  "<INFO>\t",timer(),"\t$var_file\n";
}
print  "\n<INFO>\t",timer(),"\tStart creating joint variant table logic...\n";
tbjoin($VAR_dir,$POS_OUT,$CALL_OUT,$JOIN_OUT,$group_name,$ref,$refg,$micovf,$micovr,$miphred20,$mifreq,@var_files_new);
print  "<INFO>\t",timer(),"\tFinished creating joint variant table logic!\n";
@var_files         	=	();
@var_files_new		=	();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBamend:
if($step eq 'TBamend') {
    	print  "\n<INFO>\t",timer(),"\t### [TBamend] workflow selected ###\n";
}
opendir(JOINDIR,"$JOIN_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $JOIN_OUT: $!";
opendir(AMENDDIR,"$AMEND_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $AMEND_OUT: $!";
@join_files		=	grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode3_samples$sample_number\.tab$/ && -f "$JOIN_OUT/$_" } readdir(JOINDIR);
@amend_files            =       grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode3_samples$sample_number\_amended_u$unambigous\_phylo_w$window.tab$/ && -f "$AMEND_OUT/$_" } readdir(AMENDDIR);
if(scalar(@join_files) == 0) {
    	print  "\n<ERROR>\t",timer(),"\tNo joint variants file to amend! Check content of $JOIN_OUT!\n";
    	exit 1;
}
%check_up               =       map { (my $id = $_) =~ s/\_amended_u$unambigous\_phylo_w$window.tab//; $id => $id; } @amend_files;
for(my $i = 0;$i < scalar(@join_files);$i++) {
        my $tmp         =       $join_files[$i];
        $tmp            =~      s/\.tab//;
        if(!exists $check_up{$tmp}) {
                push(@join_files_new,$join_files[$i]);
        }
}
if(scalar(@join_files_new) == 0) {
        print  "\n<ERROR>\t",timer(),"\tA Joint variant file with $group_name and identical parameter setting is already existing! Check content of $JOIN_OUT!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tAmending joint variant table:\n";
foreach my $join_file (sort { $a cmp $b } @join_files_new) {
    	print  "<INFO>\t",timer(),"\t$join_file\n";
}
print  "\n<INFO>\t",timer(),"\tStart amending joint variant table logic...\n";
tbamend($JOIN_OUT,$AMEND_OUT,$micovf,$micovr,$miphred20,$mifreq,$unambigous,$window,$categories,$resi_list_master,$int_regions,@join_files_new);
print  "<INFO>\t",timer(),"\tFinished variant calling from position table logic!\n";
@join_files		=	();
@amend_files		=	();
@join_files_new         =       ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBstrains:
if($step eq 'TBstrains') {
    	print  "\n<INFO>\t",timer(),"\t### [TBstrains] workflow selected ###\n";
}
opendir(POSDIR,"$POS_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $POS_OUT: $!";
@pos_files          	=       grep { $_ =~ /^\w.*\.gatk_position_table\.tab$/ && -f "$POS_OUT/$_" } readdir(POSDIR);
if(scalar(@pos_files) == 0) {
    	print  "\n<ERROR>\t",timer(),"\tNo files to call strains! Check content of $POS_OUT!\n";
    	exit 1;
}
print  "\n<INFO>\t",timer(),"\tCalling strains:\n";
foreach my $pos_file (sort { $a cmp $b } @pos_files) {
        print  "<INFO>\t",timer(),"\t$pos_file\n";
}
print  "\n<INFO>\t",timer(),"\tStart strain call logic...\n";
tbstrains($VAR_dir,$POS_OUT,$STRAIN_OUT,$ref,$refg,$micovf,$micovr,$mifreq,$miphred20,$date_string,@pos_files);
print  "<INFO>\t",timer(),"\tFinished strain call logic!\n";
@pos_files		=	();
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBgroups:
if($step eq 'TBgroups') {
    	print  "\n<INFO>\t",timer(),"\t### [TBgroups] workflow selected ###\n";
}
opendir(AMENDDIR,"$AMEND_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $AMEND_OUT: $!";
opendir(GROUPDIR,"$GROUPS_OUT") || die "<ERROR>\t",timer(),"\tCan\'t open directory $GROUPS_OUT: $!";
@amend_files            =       grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode3_samples$sample_number\_amended_u$unambigous\_phylo_w$window\.tab$/ && -f "$AMEND_OUT/$_" } readdir(AMENDDIR);
@group_files		=	grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode3_samples$sample_number\_amended_u$unambigous\_phylo_w$window\.matrix$/ && -f "$GROUPS_OUT/$_" } readdir(GROUPDIR);
if(scalar(@amend_files) == 0) {
    	print  "\n<ERROR>\t",timer(),"\tNo file to call groups! Check content of $AMEND_OUT and the group name $group_name!\n";
    	exit 1;
}
%check_up               =       map { (my $id = $_) =~ s/\.matrix$//; $id => $id; } @group_files;
for(my $i = 0;$i < scalar(@amend_files);$i++) {
        my $tmp         =       $amend_files[$i];
        $tmp		=~	s/\.tab//;
	if(!exists $check_up{$tmp}) {
                push(@amend_files_new,$amend_files[$i]);
        }
}
if(scalar(@amend_files_new) == 0) {
        print  "\n<ERROR>\t",timer(),"\tA group call with group name $group_name and identical parameter setting is already existing! Check content of $GROUPS_OUT!\n";
        exit 1;
}
print  "\n<INFO>\t",timer(),"\tCalling groups from:\n";
foreach my $amend_file (sort { $a cmp $b } @amend_files_new) {
    	print  "<INFO>\t",timer(),"\t$amend_file\n";
}
print  "\n<INFO>\t",timer(),"\tStart calling groups logic...\n";
tbgroups($AMEND_OUT,$GROUPS_OUT,$distance,@amend_files_new);
print  "<INFO>\t",timer(),"\tFinished calling groups logic!\n";
@amend_files            =       ();
@group_files		=	();
@amend_files_new	=	();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



print  "\n<INFO>\t",timer(),"\t### TBseq finished!!! ###\n";
