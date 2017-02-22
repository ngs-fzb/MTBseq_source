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

package TBreads;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###														###
### Description: This package converts the native Illumina naming scheme into a more practical naming scheme.	###
###														###
### Input: 	Sample-ID-LibID_SampleNumber_LaneNumber_Dir_ID.fastq.gz						###
### Output:	Sample-ID_LibID_Machine_Run_Readlength_Dir.fastq.gz						###
###														###
###################################################################################################################

$VERSION	=	1.00;
@ISA		=	qw(Exporter);
@EXPORT		= 	qw(tbreads);


sub tbreads {
	# Switches autoflush for direct printing on.
	$|			=	1;     
	# Get parameter and input from front-end.
	my $W_dir		=	shift;
	my $machine		=	shift;
	my $run_number		=	shift;
	my $seq_length		=	shift;
	$seq_length		=	$seq_length . "bp";
	my @fastq_files		=	@_;
	my %samples_hash;
	my %lib_hash;
	# Start logic...
	foreach my $file(sort {$a cmp $b } @fastq_files) {
		my $new_file		=	$file;
		my @filename            =       split(/_/,$new_file);
		if((scalar(@filename) < 5 || scalar(@filename) > 6) && (!$new_file =~ /-lib/)) {
			print  "<ERROR>\t",timer(),"\tSkipping $new_file\. Bad file name for $new_file\. Check README.pdf for more information!\n";
                        next;
                }
		@filename		=	();
		$new_file		=~	s/(.*)-lib(.*)/$1\_lib$2/;
		$new_file		=~	s/(.*)-wdh(.*)/$1$2/;
		$new_file		=~	s/(.*)-wdh2(.*)/$1$2/;
		$new_file		=~	s/(.*)-neu(.*)/$1$2/;
		move("$W_dir/$file","$W_dir/$new_file") || die "<ERROR>\t",timer(),"\tmove failed: $!\n";
		@filename		=	split(/_/,$new_file);
		my $sampleID		=	$filename[0];
		my $libID		=	$filename[1];
		my $sampleNumber	=	$filename[2];
		my $laneNumber		=	$filename[3];
		my $dir;
		for(my $i = 0;$i < scalar(@filename);$i++) {
			if($filename[$i] =~ /R\d+/) {
				$dir	=	$filename[$i];
				$dir	=~      s/.fastq.gz//;
			}
		}
		push(@{$samples_hash{$sampleID}{$sampleNumber}{$dir}},$new_file);
		$lib_hash{$sampleID}	=	$libID;
	}
	# Do the logic for every sample within the working directory.
	foreach my $sampleID(sort { $a cmp $b} keys(%samples_hash)) {
		my @sampleNumber	=	keys(%{$samples_hash{$sampleID}});
		if(scalar(@sampleNumber) > 1) {
			print  "<WARN>\t",timer(),"\tSkipping $sampleID\. More than one dataset for sample $sampleID!\n";
			next;
		}
		# Look for multiple lane files and start merging.
		foreach my $sampleNumber(sort { $a cmp $b } @sampleNumber) {
			print  "<INFO>\t",timer(),"\tStart file processing for $sampleID...\n";
			my @forward_files	=	@{$samples_hash{$sampleID}{$sampleNumber}{R1}};
			my @reverse_files	=	@{$samples_hash{$sampleID}{$sampleNumber}{R2}};
			if(scalar(@forward_files) == 1 && scalar(@reverse_files == 1)) {
				print  "<WARN>\t",timer(),"\tSkipping $sampleID\. No multiple files found!\n";
				next;
			}
			if(scalar(@forward_files) != scalar(@reverse_files)) {
				print  "<WARN>\t",timer(),"\tSkipping $sampleID\. Inconsistent number of forward/reverse files! Maybe lane missing?\n";
				next;
			}
			my $libID		=	$lib_hash{$sampleID};
			my @output_file		=	($sampleID,$libID,$machine,$run_number,$seq_length);
			my $output_file_R1	=	join("_",@output_file);
			my $output_file_R2	=	join("_",@output_file);
			$output_file_R1         .=      "_R1.fastq.gz";
			$output_file_R2		.=	"_R2.fastq.gz";
			if(-f "$W_dir/$output_file_R1") {
				print  "<WARN>\t",timer(),"\tSkipping $sampleID\. $output_file_R1 already exist!\n";
				next;
			}
			if(-f "$W_dir/$output_file_R2") {
				print  "<WARN>\t",timer(),"\tSkipping $sampleID\. $output_file_R2 already exist!\n";
				next;
			}
			my $file_string_1	=	"";
			my $file_string_2	=	"";
			foreach my $file(sort { $a cmp $b } @forward_files) {
				$file_string_1		.=		" $W_dir/$file";
			}
			foreach my $file(sort { $a cmp $b } @reverse_files) {
				$file_string_2		.=		" $W_dir/$file";
			}
			print  "<INFO>\t",timer(),"\tMerging multiple files for $sampleID and standardize file naming...\n";
			system("zcat $file_string_1 | gzip -c > $output_file_R1");
			system("zcat $file_string_2 | gzip -c > $output_file_R2");
			print  "<INFO>\t",timer(),"\tRemoving multiple files for $sampleID...\n";
			foreach my $file(sort { $a cmp $b } @forward_files) {
				unlink("$W_dir/$file") || warn "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
			}
			foreach my $file(sort { $a cmp $b } @reverse_files) {
				unlink("$W_dir/$file") || warn "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
			}
			# Finished.
			print  "<INFO>\t",timer(),"\tFinished file processing for $sampleID!\n";
		}
	}
}


1;
