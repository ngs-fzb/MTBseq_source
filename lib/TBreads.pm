#!/usr/bin/perl

# tabstop is set to 8.

package TBreads;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION	=	1.0.0;
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
	# Start logic...
	foreach my $file (sort {$a cmp $b } @fastq_files) {
		my @file_name            =       split(/_/,$file);
		if(scalar(@file_name) < 3) {
			print $logprint "<ERROR>\t",timer(),"\tRequire at least [SampleID]_[LibID]_[Direction] field within the file name of $file! Check README.pdf for more information!\n";
                        exit;
                }
		my $sampleID		=	shift(@file_name);
		my $libID		=	shift(@file_name);
		my $dir			=	pop(@file_name);
		my $machine_f		=	shift(@file_name);
		my $run_f		=	shift(@file_name);
		if(!$dir || !$dir =~ /R\d/) {
			print $logprint "<ERROR>\t",timer(),"\tCannot find [Direction] field within filename of $file! Check README.pdf for more information!\n";
			exit;
		}
		$dir			=~	s/\.fastq.*$//;
		push(@{$samples_hash{$sampleID}{$libID}{$dir}},$file);
		$machine_f_hash{$sampleID}{$libID}		=	$machine_f	if($machine_f && $run_f);
		$run_f_hash{$sampleID}{$libID}			=	$run_f		if($machine_f && $run_f);
	}
	@fastq_files = ();
	# Start logic for every sample within the working directory.
	foreach my $sampleID (sort { $a cmp $b } keys %samples_hash) {
		my @libID		=	keys(%{$samples_hash{$sampleID}});
		if(scalar(@libID) > 1) {
			print $logprint "<INFO>\t",timer(),"\t$sampleID was sequenced with more then one library!\n";
		}
		# Look for multiple sample_lib files and start merging.
		foreach my $libID (sort { $a cmp $b } @libID) {
			print $logprint "<INFO>\t",timer(),"\tStart file processing for $sampleID\_$libID...\n";
			my @forward_files	=	@{$samples_hash{$sampleID}{$libID}{R1}} if(exists $samples_hash{$sampleID}{$libID}{R1});
			my @reverse_files	=	@{$samples_hash{$sampleID}{$libID}{R2}} if(exists $samples_hash{$sampleID}{$libID}{R2});
			my @uni_files		=	@{$samples_hash{$sampleID}{$libID}{R0}} if(exists $samples_hash{$sampleID}{$libID}{R0});
			my @output_file_R1      =       ($sampleID,$libID,$machine,$run_number);
			my @output_file_R2      =       ($sampleID,$libID,$machine,$run_number);
			my @output_file_R0      =       ($sampleID,$libID,$machine,$run_number);
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
                        	print $logprint "<INFO>\t",timer(),"\tSkipping $sampleID\_$libID. Files already exists!\n";
                                next;
                        }
			if(@forward_files && @reverse_files) {
				if(scalar(@forward_files) == 1 && scalar(@reverse_files == 1)) {
					print $logprint "<INFO>\t",timer(),"\t$sampleID\_$libID has no multiple read files!\n";
				}
				if(scalar(@forward_files) != scalar(@reverse_files)) {
					print $logprint "<INFO>\t",timer(),"\tInconsistent number of forward/reverse files for $sampleID\_$libID!\n";
				}
				foreach my $file (@forward_files) {
					if($file =~ /fastq\.gz$/) {
						print $logprint "<INFO>\t",timer(),"\t$file is already compressed!\n";
					} 
					else {
						print $logprint "<INFO>\t",timer(),"\t$file is not compressed! Start compressing $file...\n";
                                                system("gzip $file");
						print $logprint "<INFO>\t",timer(),"\tFinished compressing $file!\n";
					}
				}
				foreach my $file (@reverse_files) {
                                        if($file =~ /fastq\.gz$/) {
                                                print $logprint "<INFO>\t",timer(),"\t$file is already compressed!\n";
                                        }
                                        else {
						print $logprint "<INFO>\t",timer(),"\t$file is not compressed! Start compressing $file...\n";
                                                system("gzip $file");
						print $logprint "<INFO>\t",timer(),"\tFinished compressing $file!\n";
                                        }
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
				print $logprint "<INFO>\t",timer(),"\tStart merging files for $sampleID\_$libID...\n";
				system("zcat $file_string_1 | gzip -c > $output_file_R1");
				system("zcat $file_string_2 | gzip -c > $output_file_R2");
				print $logprint "<INFO>\t",timer(),"\tFinished merging files for $sampleID\_$libID!\n";
				print $logprint "<INFO>\t",timer(),"\tStart removing multiple files for $sampleID\_$libID...\n";
				foreach my $file (sort { $a cmp $b } @forward_files) {
					$file =~ s/\.fastq.*$//;
					unlink("$W_dir/$file.fastq.gz");
				}
				foreach my $file (sort { $a cmp $b } @reverse_files) {
					$file =~ s/\.fastq.*$//;
					unlink("$W_dir/$file.fastq.gz");
				}
				print $logprint "<INFO>\t",timer(),"\tFinished removing multiple files for $sampleID\_$libID!\n";
			}
			@forward_files = ();
			@reverse_files = ();
			if(@uni_files) {
				if(scalar(@uni_files) == 1) {
					print $logprint "<INFO>\t",timer(),"\t$sampleID\_$libID has no multiple files!\n";               
				}
				foreach my $file (@uni_files) {
					if($file =~ /fastq\.gz$/) {
                                                print $logprint "<INFO>\t",timer(),"\t$file already compressed!\n";
                                        }
                                        else {
						print $logprint "<INFO>\t",timer(),"\t$file is not compressed! Start compressing $file...\n";
                                                system("gzip $file");
						print $logprint "<INFO>\t",timer(),"\tFinished compressing $file!\n";
					}
                                } 
				my $file_string_1       	=       "";
				foreach my $file (sort { $a cmp $b } @uni_files) {
					$file =~ s/\.fastq.*$//;					
                                        $file_string_1          .=	" $W_dir/$file.fastq.gz";
                                }
				print $logprint "<INFO>\t",timer(),"\tStart merging files for $sampleID\_$libID...\n";
				system("zcat $file_string_1 | gzip -c > $output_file_R0");
				print $logprint "<INFO>\t",timer(),"\tFinished merging files for $sampleID\_$libID!\n";
				print $logprint "<INFO>\t",timer(),"\tRemoving multiple files for $sampleID\_$libID...\n";
                                foreach my $file (sort { $a cmp $b } @uni_files) {
					$file =~ s/\.fastq.*$//;
                                        unlink("$W_dir/$file.fastq.gz");
                                }
				print $logprint "<INFO>\t",timer(),"\tFinished removing multiple files for $sampleID\_$libID!\n";
			}
			@uni_files = ();
		}
		print $logprint "<INFO>\t",timer(),"\tFinished read file processing for $sampleID!\n";
		@libID = ();
	}
	undef(%samples_hash);
	undef(%machine_f_hash);
	undef(%run_f_hash);
}


1;
