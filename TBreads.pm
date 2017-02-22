#!/usr/bin/perl

=head1

        TBseq - a computational pipeline for detecting variants in NGS-data

        Copyright (C) 2016 Thomas A. Kohl, Robin Koch, Maria R. De Filippo, Viola Schleusener, Christian Utpatel, Daniela M. Cirillo, Stefan Niemann

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

# tabstop is set to 8.

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

$VERSION	=	1.10;
@ISA		=	qw(Exporter);
@EXPORT		= 	qw(tbreads);


sub tbreads {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $W_dir		=	shift;
	my $machine		=	shift;
	my $run_number		=	shift;
	my $seq_length		=	shift;
	$seq_length		=	$seq_length;
	my @fastq_files		=	@_;
	my %samples_hash;
	# Start logic...
	foreach my $file(sort {$a cmp $b } @fastq_files) {
		my $new_file		=	$file;
		my $new_file2		=	$file;
		my @filename            =       split(/_/,$new_file);
		if((scalar(@filename) < 3 || scalar(@filename) > 6) && (!$new_file =~ /lib/)) {
			print $logprint "<WARN>\t",timer(),"\tSkipping $file! Require at least [SampleID]_[LibID]_[Direction] field wihtin the file name! Check README.pdf for more information!\n";
                        next;
                }
		$new_file		=~	s/(.*)-lib(.*)/$1\_lib$2/;
		$new_file		=~	s/(.*)-wdh(.*)/$1$2/;
		$new_file		=~	s/(.*)-wdh2(.*)/$1$2/;
		$new_file		=~	s/(.*)-neu(.*)/$1$2/;
		move("$W_dir/$file","$W_dir/$new_file") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
		@filename		=	split(/_/,$new_file);
		my $sampleID		=	shift(@filename);
		my $libID		=	shift(@filename);
		my $dir;
		for(my $i = 0;$i < scalar(@filename);$i++) {
			if($filename[$i] =~ /R\d+/) {
				$dir	=	$filename[$i];
			}
		}
		if(!$dir) {
			print $logprint "<WARN>\t",timer(),"\tSkipping $new_file! Cannot find [Direction] field in file name! Check README.pdf for more information!\n";
			next;
		}
		if(scalar(@filename) == 1) {
			$new_file2	=	"$sampleID\_$libID\_$machine\_$run_number\_$seq_length\_$dir";
			move("$W_dir/$new_file","$W_dir/$new_file2") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
			$dir                    =~      s/\.fastq.*$//;
			push(@{$samples_hash{$sampleID}{$libID}{$dir}},$new_file2);
		}
		else {
			$dir			=~	s/\.fastq.*$//;
			push(@{$samples_hash{$sampleID}{$libID}{$dir}},$new_file);
		}
	}
	# Do the logic for every sample within the working directory.
	foreach my $sampleID(sort { $a cmp $b } keys(%samples_hash)) {
		my @libID	=	keys(%{$samples_hash{$sampleID}});
		if(scalar(@libID) > 1) {
			print $logprint "<INFO>\t",timer(),"\t$sampleID was sequenced with more then one library!\n";
		}
		# Look for multiple lane files and start merging.
		foreach my $libID(sort { $a cmp $b } @libID) {
			print $logprint "<INFO>\t",timer(),"\tStart file processing for $sampleID...\n";
			my @forward_files	=	@{$samples_hash{$sampleID}{$libID}{R1}} if(exists $samples_hash{$sampleID}{$libID}{R1});
			my @reverse_files	=	@{$samples_hash{$sampleID}{$libID}{R2}} if(exists $samples_hash{$sampleID}{$libID}{R2});
			my @uni_files		=	@{$samples_hash{$sampleID}{$libID}{R0}} if(exists $samples_hash{$sampleID}{$libID}{R0});
			if(@forward_files && @reverse_files) {
				if(scalar(@forward_files) == 1 && scalar(@reverse_files == 1)) {
					print $logprint "<INFO>\t",timer(),"\t$sampleID has no multiple files but will be compressed anyway!\n";
					if($forward_files[0] =~ /fastq.gz$/) {
						next;
					}
					else {
						system("gzip $forward_files[0]");
						print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $forward_files[0]...\n";
						unlink("$W_dir/$forward_files[0]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $forward_files[0]: No such file!\n";
						print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $forward_files[0]!\n"; 
						next;
					}
					if($reverse_files[0] =~ /fastq.gz$/) {
						next;
					}
					else {
						system("gzip $reverse_files[0]");
						print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $reverse_files[0]...\n";
						unlink("$W_dir/$reverse_files[0]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $reverse_files[0]: No such file!\n";
						print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $reverse_files[0]!\n";
						next;
					}
					next;
				}
				if(scalar(@forward_files) != scalar(@reverse_files)) {
					print $logprint "<WARN>\t",timer(),"\tSkipping $sampleID! Inconsistent number of forward/reverse files! Maybe a lane is missing for one direction?\n";
					next;
				}
				my @output_file		=	($sampleID,$libID,$machine,$run_number,$seq_length);
				my $output_file_R1	=	join("_",@output_file);
				my $output_file_R2	=	join("_",@output_file);
				$output_file_R1         .=      "_R1.fastq.gz";
				$output_file_R2		.=	"_R2.fastq.gz";
				if(-f "$W_dir/$output_file_R1") {
					print $logprint "<WARN>\t",timer(),"\tSkipping $sampleID, $output_file_R1 already exists!\n";
					next;
				}
				if(-f "$W_dir/$output_file_R2") {
					print $logprint "<WARN>\t",timer(),"\tSkipping $sampleID, $output_file_R2 already exists!\n";
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
				print $logprint "<INFO>\t",timer(),"\tMerging and/or compressing files for $sampleID and standardize file naming...\n";
				system("zcat $file_string_1 | gzip -c > $output_file_R1");
				system("zcat $file_string_2 | gzip -c > $output_file_R2");
				print $logprint "<INFO>\t",timer(),"\tRemoving multiple files for $sampleID...\n";
				foreach my $file(sort { $a cmp $b } @forward_files) {
					unlink("$W_dir/$file") || print $logprint "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
				}
				foreach my $file(sort { $a cmp $b } @reverse_files) {
					unlink("$W_dir/$file") || print $logprint "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
				}
			}
			if(@uni_files) {
				if(scalar(@uni_files) == 1) {
					print $logprint "<INFO>\t",timer(),"\t$sampleID has no multiple files but will be compressed anyway!\n";
                                        if($uni_files[0] =~ /fastq.gz$/) {
                                        	next;
					}
                                        else {
                                                system("gzip $uni_files[0]");
					        print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $uni_files[0]...\n";
                                                unlink("$W_dir/$uni_files[0]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $uni_files[0]: No such file!\n";
                                                print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $uni_files[0]!\n";
						next;					
					}
				}
				my @output_file         =       ($sampleID,$libID,$machine,$run_number,$seq_length);
				my $output_file_R0      =       join("_",@output_file);
				$output_file_R0         .=      "_R0.fastq.gz";
                                if(-f "$W_dir/$output_file_R0") {
                                        print $logprint "<WARN>\t",timer(),"\tSkipping $sampleID, $output_file_R0 already exists!\n";
                                        next;
                                }
                                my $file_string_1       =       "";
				foreach my $file(sort { $a cmp $b } @uni_files) {
                                        $file_string_1          .=              " $W_dir/$file";
                                }
				print $logprint "<INFO>\t",timer(),"\tMerging and/or compressing files for $sampleID and standardize file naming...\n";
				system("zcat $file_string_1 | gzip -c > $output_file_R0");
				print $logprint "<INFO>\t",timer(),"\tRemoving multiple files for $sampleID...\n";
                                foreach my $file(sort { $a cmp $b } @uni_files) {
                                        unlink("$W_dir/$file") || print $logprint "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
                                }
			}
			print $logprint "<INFO>\t",timer(),"\tFinished file processing for $sampleID!\n";
		}
	}
}


1;
