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
### Input: 	[SampleID]_[LibID]_[Dir].fastq.gz								###
### Output:	[SampleID]_[LibID]_[Machine]_[Run]_[Dir].fastq.gz						###
###														###
###################################################################################################################

$VERSION	=	1.15;
@ISA		=	qw(Exporter);
@EXPORT		= 	qw(tbreads);


sub tbreads {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $W_dir		=	shift;
	my $machine		=	shift;
	my $run_number		=	shift;
	my @fastq_files		=	@_;
	my %samples_hash;
	my %machine_f_hash;
	my %run_f_hash;
	my %seq_f_hash;
	# Start logic...
	foreach my $file (sort {$a cmp $b } @fastq_files) {
		my @file_name            =       split(/_/,$file);
		if(scalar(@file_name) < 3 || (!$file =~ /lib/)) {
			print $logprint "<WARN>\t",timer(),"\tRequire at least [SampleID]_[LibID]_[Direction] field within the file name of $file! [LibID] must start with \'lib\' Check README.pdf for more information!\n";
                        exit;
                }
		my $sampleID		=	shift(@file_name);
		my $libID		=	shift(@file_name);
		my $dir			=	pop(@file_name);
		my $machine_f		=	shift(@file_name);
		my $run_f		=	shift(@file_name);
		my $seq_f               =       shift(@file_name);
		if(!$dir || !$dir =~ /R\d/) {
			print $logprint "<WARN>\t",timer(),"\tCannot find [Direction] field in file name of $file! Check README.pdf for more information!\n";
			exit;
		}
		$dir			=~	s/\.fastq.*$//;
		push(@{$samples_hash{$sampleID}{$libID}{$dir}},$file);
		$machine_f_hash{$sampleID}{$libID}		=	$machine_f	if($machine_f && $run_f && $seq_f);
		$run_f_hash{$sampleID}{$libID}			=	$run_f		if($machine_f && $run_f && $seq_f);
		$seq_f_hash{$sampleID}{$libID}{$dir}          	=	$seq_f		if($machine_f && $run_f && $seq_f);
	}
	@fastq_files = ();
	# Do the logic for every sample within the working directory.
	foreach my $sampleID (sort { $a cmp $b } keys %samples_hash) {
		my @libID		=	keys(%{$samples_hash{$sampleID}});
		if(scalar(@libID) > 1) {
			print $logprint "<INFO>\t",timer(),"\t$sampleID was sequenced with more then one library!\n";
		}
		# Look for multiple sample_lib files and start merging.
		foreach my $libID (sort { $a cmp $b } @libID) {
			print $logprint "<INFO>\t",timer(),"\tStart file processing for $sampleID...\n";
			my @forward_files	=	@{$samples_hash{$sampleID}{$libID}{R1}} if(exists $samples_hash{$sampleID}{$libID}{R1});
			my @reverse_files	=	@{$samples_hash{$sampleID}{$libID}{R2}} if(exists $samples_hash{$sampleID}{$libID}{R2});
			my @uni_files		=	@{$samples_hash{$sampleID}{$libID}{R0}} if(exists $samples_hash{$sampleID}{$libID}{R0});
			my $seq_len_R1		=	"1bp";
			my $seq_len_R2          =       "1bp";
			my $seq_len_R0          =       "1bp";
			$seq_len_R1		=	$seq_f_hash{$sampleID}{$libID}{R1}	if(exists $seq_f_hash{$sampleID}{$libID}{R1});
			$seq_len_R2		=	$seq_f_hash{$sampleID}{$libID}{R2}      if(exists $seq_f_hash{$sampleID}{$libID}{R2});
                        $seq_len_R0		=       $seq_f_hash{$sampleID}{$libID}{R0}      if(exists $seq_f_hash{$sampleID}{$libID}{R0});
			my @output_file_R1      =       ($sampleID,$libID,$machine,$run_number,$seq_len_R1);
			my @output_file_R2      =       ($sampleID,$libID,$machine,$run_number,$seq_len_R2);
			my @output_file_R0      =       ($sampleID,$libID,$machine,$run_number,$seq_len_R0);
			$output_file_R1[2]	=	$machine_f_hash{$sampleID}{$libID}	if(exists $machine_f_hash{$sampleID}{$libID});
			$output_file_R2[2]      =       $machine_f_hash{$sampleID}{$libID}      if(exists $machine_f_hash{$sampleID}{$libID});
			$output_file_R0[2]      =       $machine_f_hash{$sampleID}{$libID}      if(exists $machine_f_hash{$sampleID}{$libID});
                        $output_file_R1[3]      =       $run_f_hash{$sampleID}{$libID}      	if(exists $run_f_hash{$sampleID}{$libID});
                        $output_file_R2[3]      =       $run_f_hash{$sampleID}{$libID}      	if(exists $run_f_hash{$sampleID}{$libID});
                        $output_file_R0[3]      =       $run_f_hash{$sampleID}{$libID}      	if(exists $run_f_hash{$sampleID}{$libID});                        
			my $output_file_R1      =       join("_",@output_file_R1);
                        my $output_file_R2      =       join("_",@output_file_R2);
                        my $output_file_R0      =       join("_",@output_file_R0);
			$output_file_R1         .=      "_R1.fastq.gz";
                        $output_file_R2         .=      "_R2.fastq.gz";
			$output_file_R0         .=      "_R0.fastq.gz";
                        if((-f "$W_dir/$output_file_R1" && -f "$W_dir/$output_file_R2") || -f "$W_dir/$output_file_R0") {
                        	print $logprint "<INFO>\t",timer(),"\tSkipping $sampleID, files already exists!\n";
                                next;
                        }
			if(@forward_files && @reverse_files) {
				my $seq_len	=	get_seq_len($logprint,$W_dir,@forward_files);
				$output_file_R1 =~ 	s/(.*)_.*bp.*$/$1\_$seq_len\_R1.fastq.gz/;
				$seq_len        =       get_seq_len($logprint,$W_dir,@reverse_files);
				$output_file_R2 =~	s/(.*)_.*bp.*$/$1\_$seq_len\_R2.fastq.gz/;
				if(scalar(@forward_files) == 1 && scalar(@reverse_files == 1)) {
					print $logprint "<INFO>\t",timer(),"\t$sampleID has no multiple files!\n";
					if($forward_files[0] =~ /fastq\.gz$/) {
						print $logprint "<INFO>\t",timer(),"\t$forward_files[0] already compressed! Will only rename the file.\n";
						move("$W_dir/$forward_files[0]","$output_file_R1") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
					} 
					else {
						print $logprint "<INFO>\t",timer(),"\t$forward_files[0] is not compressed! Will compress and rename the file!\n";
                                                system("gzip $forward_files[0]");
                                                print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $forward_files[0]...\n";
                                                unlink("$W_dir/$forward_files[0]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $forward_files[0]: No such file!\n";
                                                print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $forward_files[0]!\n";
						move("$W_dir/$forward_files[0].gz","$output_file_R1") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
					}
                                        if($reverse_files[0] =~ /fastq\.gz$/) {
                                                print $logprint "<INFO>\t",timer(),"\t$reverse_files[0] already compressed! Will only rename the file.\n";
						move("$W_dir/$reverse_files[0]","$output_file_R2") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
                                        }
                                        else {
						print $logprint "<INFO>\t",timer(),"\t$reverse_files[0] is not compressed! Will compress and rename the file!\n";
                                                system("gzip $reverse_files[0]");
                                                print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $reverse_files[0]...\n";
                                                unlink("$W_dir/$reverse_files[0]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $reverse_files[0]: No such file!\n";
                                                print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $reverse_files[0]!\n";
						move("$W_dir/$reverse_files[0].gz","$output_file_R2") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
                                        }
				}
				else {
					if(scalar(@forward_files) != scalar(@reverse_files)) {
						print $logprint "<INFO>\t",timer(),"\tInconsistent number of forward/reverse files for $sampleID!\n";
                                	}
					for(my $i = 0; $i < scalar(@forward_files); $i++) {
						if($forward_files[$i] =~ /fastq\.gz$/) {
							print $logprint "<INFO>\t",timer(),"\t$forward_files[$i] already compressed!\n";
						}
						else {
							print $logprint "<INFO>\t",timer(),"\t$forward_files[$i] is not compressed! Will compress and rename the file!\n";
							system("gzip $forward_files[$i]");
							print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $forward_files[$i]...\n";
							unlink("$W_dir/$forward_files[$i]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $forward_files[$i]: No such file!\n";
							print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $forward_files[$i]!\n"; 
						}
					}
					for(my $i = 0; $i < scalar(@reverse_files); $i++) {
						if($reverse_files[$i] =~ /fastq\.gz$/) {
							print $logprint "<INFO>\t",timer(),"\t$reverse_files[$i] already compressed!\n";
						}
						else {
							print $logprint "<INFO>\t",timer(),"\t$reverse_files[$i] is not compressed! Will compress and rename the file!\n";
							system("gzip $reverse_files[$i]");
							print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $reverse_files[$i]...\n";
							unlink("$W_dir/$reverse_files[$i]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $reverse_files[$i]: No such file!\n";
							print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $reverse_files[$i]!\n";
						}
					}
				}
                        	if((-f "$W_dir/$output_file_R1" && -f "$W_dir/$output_file_R2") || -f "$W_dir/$output_file_R0") {
                               		next;
                        	}				
				my $file_string_1	=	"";
				my $file_string_2	=	"";
				foreach my $file (sort { $a cmp $b } @forward_files) {
					$file =~ s/\.fastq.*$//;
					$file_string_1	.=	" $W_dir/$file.fastq.gz";
				}
				foreach my $file (sort { $a cmp $b } @reverse_files) {
					$file =~ s/\.fastq.*$//;
					$file_string_2	.=	" $W_dir/$file.fastq.gz";
				}
				print $logprint "<INFO>\t",timer(),"\tStart merging and compressing files for $sampleID...\n";
				system("zcat $file_string_1 | gzip -c > $output_file_R1");
				system("zcat $file_string_2 | gzip -c > $output_file_R2");
				print $logprint "<INFO>\t",timer(),"\tFinished merging and compressing files for $sampleID!\n";
				print $logprint "<INFO>\t",timer(),"\tStart removing multiple files for $sampleID...\n";
				foreach my $file (sort { $a cmp $b } @forward_files) {
					$file =~ s/\.fastq.*$//;
					unlink("$W_dir/$file.fastq.gz") || print $logprint "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
				}
				foreach my $file (sort { $a cmp $b } @reverse_files) {
					$file =~ s/\.fastq.*$//;
					unlink("$W_dir/$file.fastq.gz") || print $logprint "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
				}
				print $logprint "<INFO>\t",timer(),"\tFinished removing multiple files for $sampleID!\n";
			}
			@forward_files = ();
			@reverse_files = ();
			if(@uni_files) {
				my $seq_len 	=	get_seq_len($logprint,$W_dir,@uni_files);
				$output_file_R0 =~ 	s/(.*)_.*bp.*$/$1\_$seq_len\_R0.fastq.gz/;
				if(scalar(@uni_files) == 1) {
					print $logprint "<INFO>\t",timer(),"\t$sampleID has no multiple files!\n";
                                        if($uni_files[0] =~ /fastq\.gz$/) {
                                                print $logprint "<INFO>\t",timer(),"\t$uni_files[0] already compressed! Will only rename the file.\n";
                                                move("$W_dir/$uni_files[0]","$output_file_R0") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
                                        }
                                        else {
						print $logprint "<INFO>\t",timer(),"\t$uni_files[0] is not compressed! Will compress and rename the file!\n";
                                                system("gzip $uni_files[0]");
                                                print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $uni_files[0]...\n";
                                                unlink("$W_dir/$uni_files[0]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $uni_files[0]: No such file!\n";
                                                print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $uni_files[0]!\n";
                                                move("$W_dir/$uni_files[0].gz","$output_file_R0") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
                                        }
                                } 
				else {
					for(my $i = 0; $i < scalar(@uni_files); $i++) {
						if($uni_files[$i] =~ /fastq.gz$/) {
                	         			print $logprint "<INFO>\t",timer(),"\t$uni_files[$i] already compressed!\n";               	
						}
						else {
							print $logprint "<INFO>\t",timer(),"\t$uni_files[$i] is not compressed! Will compress and rename the file!\n";
                                                	system("gzip $uni_files[$i]");
					        	print $logprint "<INFO>\t",timer(),"\tStart removing uncompressed $uni_files[$i]...\n";
                                                	unlink("$W_dir/$uni_files[$i]") || print $logprint "<WARN>\t",timer(),"\tCan't delete $uni_files[$i]: No such file!\n";
                                                	print $logprint "<INFO>\t",timer(),"\tFinished removing uncompressed $uni_files[$i]!\n";
						}
					}
				}
                        	if((-f "$W_dir/$output_file_R1" && -f "$W_dir/$output_file_R2") || -f "$W_dir/$output_file_R0") {
                                	next;
                        	}                                
				my $file_string_1       	=       "";
				foreach my $file (sort { $a cmp $b } @uni_files) {
					$file =~ s/\.fastq.*$//;					
                                        $file_string_1          .=	" $W_dir/$file.fastq.gz";
                                }
				print $logprint "<INFO>\t",timer(),"\tStart merging and compressing files for $sampleID...\n";
				system("zcat $file_string_1 | gzip -c > $output_file_R0");
				print $logprint "<INFO>\t",timer(),"\tRemoving multiple files for $sampleID...\n";
				print $logprint "<INFO>\t",timer(),"\tFinished merging and compressing files for $sampleID!\n";
                                foreach my $file (sort { $a cmp $b } @uni_files) {
					$file =~ s/\.fastq.*$//;
                                        unlink("$W_dir/$file.fastq.gz") || print $logprint "<WARN>\t",timer(),"\tCan't delete $file: No such file!\n";
                                }
			}
			@uni_files = ();
			print $logprint "<INFO>\t",timer(),"\tFinished file processing for $sampleID!\n";
		}
		@libID = ();
	}
	undef(%samples_hash);
	undef(%machine_f_hash);
	undef(%run_f_hash);
	undef(%seq_f_hash);
}


1;
