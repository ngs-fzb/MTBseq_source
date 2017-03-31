#!/usr/bin/perl

# tabstop is set to 8.

package TBmerge;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use List::Util qw(max);
use vars qw($VERSION @ISA @EXPORT);

$VERSION	=	1.0.0;
@ISA		=	qw(Exporter);
@EXPORT		=	qw(tbmerge);

sub tbmerge {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $SAMBAMBA_dir	=	shift;
	my $BAM_OUT		=	shift;
    	my $MBAM_OUT		=	shift;
    	my $GATK_OUT		=	shift;
    	my $MPILE_OUT		=	shift;
    	my $POS_OUT		=	shift;
    	my $CALL_OUT		=	shift;
    	my $threads		=	shift;
	my @bam_files		=	@_;
	my %input;
	# Prepare input.
	foreach my $file (sort { $a cmp $b } @bam_files) {
		my @file_name		=	split(/_/,$file);
		my $sampleID		=	shift(@file_name);
		my $libID		=	shift(@file_name);
		my $machine		=	shift(@file_name);
		my $run			=	shift(@file_name);
		my $samplelib		=	$sampleID."_".$libID;
		push(@{$input{$samplelib}},$file);
	}
	@bam_files = ();
	# Start logic...
	foreach my $samplelib (sort { $a cmp $b } keys %input) {
		my @bams		=	@{$input{$samplelib}};
		if(scalar(@bams) == 1) {
			print $logprint  "<INFO>\t",timer(),"\tSkipping merging workflow for $samplelib, nothing to merge!\n";
			next;
		}
		my $sampleID			=	"";
		my $libID			=	"";
		my $machine			=	"";
		my $run				=	"";
		my $multi_merge 		= 	0;
		my $old_multi_mergelog_file	=	"";
		foreach my $file (sort { $a cmp $b } @bams) {
			my @file_name		=	split(/_/,$file);
			$sampleID		= 	shift(@file_name);
			$libID			= 	shift(@file_name);
			$machine		= 	shift(@file_name);
			$run			= 	shift(@file_name);
			if($machine =~ /^multi(\d+)/) {	
				$multi_merge			=	$1;
				$old_multi_mergelog_file	= 	$file;
				$old_multi_mergelog_file	=~	s/\.bam$/\.mergelog/;
			}
		}
		my $multisource					=	"multi" . scalar(@bams);
		if($multi_merge > 0) {
			my $correct_multisource_number		=	scalar(@bams) -1 + $multi_merge;
			$multisource				= 	"multi" . $correct_multisource_number;
		}
		my $date_string			=	timer();
		$date_string			=~	s/\[//;
		$date_string                    =~      s/\]//;
		$date_string                    =~      s/\-/_/g;
		$date_string                    =~      s/\:/_/g;
		$date_string                    =~      s/\s.*//;
		$date_string			=~	s/_/-/g;
		my $multi_file			=	$sampleID . "_" . $libID . "_" . $multisource . "_" . $date_string . ".bam";
		my $logfile			=	$sampleID . "_" . $libID . "_" . $multisource . "_" . $date_string . ".mergelog";
		unlink("$BAM_OUT/$logfile");
		if(-f "$BAM_OUT/$old_multi_mergelog_file") {
			cat($logprint,"$BAM_OUT/$old_multi_mergelog_file","$BAM_OUT/$logfile")			|| die print $logprint "<ERROR>\t",timer(),"\tcat failed: TBmerge line: ", __LINE__ , " \n";
			move("$BAM_OUT/$old_multi_mergelog_file","$MBAM_OUT/$old_multi_mergelog_file") 		|| die print $logprint "<ERROR>\t",timer(),"\tmove failed: TBmerge line: " , __LINE__ , " \n";
		}
		my $bam_string;
		foreach my $bam (sort { $a cmp $b } @bams) {
			$bam_string		.=	" $BAM_OUT/$bam";
		}
		# Merge .bam files using combined header.
		print $logprint "<INFO>\t",timer(),"\tStart merging using sambamba...\n";
		print $logprint "<INFO>\t",timer(),"\t$SAMBAMBA_dir/sambamba merge -t $threads $BAM_OUT/$multi_file $bam_string\n";
		system("$SAMBAMBA_dir/sambamba merge -t $threads $BAM_OUT/$multi_file $bam_string");
		print $logprint "<INFO>\t",timer(),"\tFinished merging using sambamba!\n";
		# Recreate index with sambamba.
		print $logprint "<INFO>\t",timer(),"\tStart recreating index using sambamba...\n";
		print $logprint "<INFO>\t",timer(),"\t$SAMBAMBA_dir/sambamba index $BAM_OUT/$multi_file\n";
		system("$SAMBAMBA_dir/sambamba index $BAM_OUT/$multi_file");
		print $logprint "<INFO>\t",timer(),"\tFinished recreating index using sambamba!\n";
		# Moving into specific directories.
		print $logprint "<INFO>\t",timer(),"\tMoving files to certain directories...\n";
		foreach my $file (sort { $a cmp $b } @bams) {
			my $bamlog_file		=	$file . "log";
			my $bai_file		=	$file . ".bai";
			cat($logprint,"$BAM_OUT/$bamlog_file","$BAM_OUT/$logfile")	|| die print $logprint "<ERROR>\t",timer(),"\tcat failed: TBmerge line: ", __LINE__ , " \n";
			move("$BAM_OUT/$file","$MBAM_OUT/$file")			|| die print $logprint "<ERROR>\t",timer(),"\tmove failed: TBmerge line: ", __LINE__ , " \n";
			move("$BAM_OUT/$bamlog_file","$MBAM_OUT/$bamlog_file")		|| die print $logprint "<ERROR>\t",timer(),"\tmove failed: TBmerge line: ", __LINE__ , " \n";
			move("$BAM_OUT/$bai_file","$MBAM_OUT/$bai_file")		|| die print $logprint "<ERROR>\t",timer(),"\tmove failed: TBmerge line: ", __LINE__ , " \n";
		}
		# Removing existing downstream files.
		print $logprint "<INFO>\t",timer(),"\tRemoving existing downstream files...\n";
		unlink glob("$GATK_OUT/$samplelib*");
		unlink glob("$MPILE_OUT/$samplelib*");
		unlink glob("$POS_OUT/$samplelib*");
		unlink glob("$CALL_OUT/$samplelib*");
		# Finished.
        	print $logprint "<INFO>\t",timer(),"\tMerging finished for $samplelib!\n";
		@bams = ();
	}
	undef(%input);
}


1;
