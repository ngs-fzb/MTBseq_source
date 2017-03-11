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

package TBbwa;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###														###
### Description: This package use BWA for mapping the reads to the reference, using SAMTOOLS for duplicate read ###
### removal and indexing of the resulting mapping file.								###
###														###
### Input:	.fastq.gz											###
### Output:	.bam, .bam.bai, .bamlog										###
###														###
###################################################################################################################

$VERSION	= 	1.12;
@ISA 		= 	qw(Exporter);
@EXPORT 	= 	qw(tbbwa);


sub tbbwa {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $W_dir               =       shift;
        my $VAR_dir             =       shift;
	my $BWA_dir             =       shift;
        my $SAMTOOLS_dir        =       shift;
	my $BAM_OUT             =       shift;
	my $ref			=	shift;
    	my $threads		=	shift;
	my @fastq_files		= 	@_;
	my %input;
	# Start logic...
	foreach my $file (sort { $a cmp $b } @fastq_files) {
		my @file_name		= 	split(/_/,$file);
		my $sampleID		= 	shift(@file_name);
		my $libID		= 	shift(@file_name);
		my $machine		= 	shift(@file_name);
		my $run			=	shift(@file_name);
		my $seqlength		=	shift(@file_name);
		my $dir			=	shift(@file_name);
		$dir			=~	s/\.fastq.gz$//;
		my $fullID		=	join("_",($sampleID,$libID,$machine,$run,$seqlength));
		$input{$fullID}{$dir}{fastq}	=	$file;
		$input{$fullID}{$dir}{sampleID}	= 	$sampleID;
		$input{$fullID}{$dir}{libID}	= 	$libID;
		$input{$fullID}{$dir}{machine}	= 	$machine;
		$input{$fullID}{$dir}{run}	= 	$run;
		$input{$fullID}{$dir}{seqlength}      =       $seqlength;
	}
	@fastq_files = ();
	foreach my $fullID (sort { $a cmp $b } keys %input) {
		my $sampleID;
		my $libID;
		my $files_string	=       "";
		my @dirs		=	sort(keys %{$input{$fullID}});
		if(scalar(@dirs) > 2) {
			print $logprint "<WARN>\t",timer(),"\tSkipping $fullID, more than two files for $fullID!\n";
			next;
		}
		foreach my $dir (sort { $a cmp $b } @dirs) {
			my $file		=	$input{$fullID}{$dir}{fastq};
			$sampleID		=	$input{$fullID}{$dir}{sampleID};
			$libID			=	$input{$fullID}{$dir}{libID};
			$files_string		.=	" $W_dir/$file";
		}
		@dirs 			=	();
		my $read_naming_scheme  =       "\'\@RG\\tID:$fullID\\tSM:$sampleID\\tPL:Illumina\\tLB:$libID\'";
		my $logfile             =       $fullID . ".bamlog";
		unlink("$BAM_OUT/$logfile") || print $logprint "<WARN>\t",timer(),"\tCan't delete $logfile: No such file!\n";;
		print $logprint "<INFO>\t",timer(),"\tFound two files for $fullID!\n";
		# Index reference
		print $logprint "<INFO>\t",timer(),"\tStart indexing reference genome $ref...\n";
		print $logprint "<INFO>\t",timer(),"\t$BWA_dir/bwa index $VAR_dir/$ref >> $BAM_OUT/$logfile\n";
		system("$BWA_dir/bwa index $VAR_dir/$ref 2>> $BAM_OUT/$logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished indexing reference genome $ref!\n";
		# Map reads with bwa-mem and -t parameter.
		print $logprint  "<INFO>\t",timer(),"\tStart BWA mapping for $fullID...\n";
		print $logprint "<INFO>\t",timer(),"\t$BWA_dir/bwa mem -t $threads -R $read_naming_scheme $VAR_dir/$ref $files_string > $BAM_OUT/$fullID.sam 2>> $BAM_OUT/$logfile\n";
		system("$BWA_dir/bwa mem -t $threads -R $read_naming_scheme $VAR_dir/$ref $files_string > $BAM_OUT/$fullID.sam 2>> $BAM_OUT/$logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished BWA mapping for $fullID!\n";
		# Convert from .sam to .bam format with -S (sam input) and (-b bam output) and -T (reference).
		print $logprint "<INFO>\t",timer(),"\tStart using samtools to convert from .sam to .bam for $fullID...\n";
		print $logprint "<INFO>\t",timer(),"\t$SAMTOOLS_dir/samtools view -@ $threads -b -T $VAR_dir/$ref -o $BAM_OUT/$fullID.bam $BAM_OUT/$fullID.sam 2>> $BAM_OUT/$logfile\n";
		system("$SAMTOOLS_dir/samtools view -@ $threads -b -T $VAR_dir/$ref -o $BAM_OUT/$fullID.bam $BAM_OUT/$fullID.sam 2>> $BAM_OUT/$logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished file conversion for $fullID!\n";
		# Sort with samtools.
		print $logprint "<INFO>\t",timer(),"\tStart using samtools for sorting of $fullID...\n";
		print $logprint "<INFO>\t",timer(),"\t$SAMTOOLS_dir/samtools sort -@ $threads -T /tmp/$fullID.sorted -o $BAM_OUT/$fullID.sorted.bam $BAM_OUT/$fullID.bam 2>> $BAM_OUT/$logfile\n";
		system("$SAMTOOLS_dir/samtools sort -@ $threads -T /tmp/$fullID.sorted -o $BAM_OUT/$fullID.sorted.bam $BAM_OUT/$fullID.bam 2>> $BAM_OUT/$logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished using samtools for sorting of $fullID!\n";
		# Indexing with samtools.
		print $logprint "<INFO>\t",timer(),"\tStart using samtools for indexing of $fullID...\n";
		print $logprint "<INFO>\t",timer(),"\t$SAMTOOLS_dir/samtools index -b $BAM_OUT/$fullID.sorted.bam 2>> $BAM_OUT/$logfile\n";
		system("$SAMTOOLS_dir/samtools index -b $BAM_OUT/$fullID.sorted.bam 2>> $BAM_OUT/$logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished using samtools for indexing of $fullID!\n";
		# Removing duplicates.
		print $logprint "<INFO>\t",timer(),"\tStart removing putative PCR duplicates from $fullID...\n";
		print $logprint "<INFO>\t",timer(),"\t$SAMTOOLS_dir/samtools rmdup $BAM_OUT/$fullID.sorted.bam $BAM_OUT/$fullID.nodup.bam 2>> $BAM_OUT/$logfile\n";
		system("$SAMTOOLS_dir/samtools rmdup $BAM_OUT/$fullID.sorted.bam $BAM_OUT/$fullID.nodup.bam 2>> $BAM_OUT/$logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished removing putative PCR duplicates for $fullID!\n";
		# Recreate index.
		print $logprint "<INFO>\t",timer(),"\tStart recreating index for $fullID...\n";
		print $logprint "<INFO>\t",timer(),"\t$SAMTOOLS_dir/samtools index -b $BAM_OUT/$fullID.nodup.bam 2>> $BAM_OUT/$logfile\n";
		system("$SAMTOOLS_dir/samtools index -b $BAM_OUT/$fullID.nodup.bam 2>> $BAM_OUT/$logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished recreating index for $fullID!\n";
		# Removing temporary files.
		print $logprint "<INFO>\t",timer(),"\tRemoving temporary files...\n";
		unlink("$BAM_OUT/$fullID.sam")			|| print $logprint "<WARN>\t",timer(),"\tCan't delete $fullID\.sam: No such file!\n";
		unlink("$BAM_OUT/$fullID.bam") 			|| print $logprint "<WARN>\t",timer(),"\tCan't delete $fullID\.bam: No such file!\n";
		unlink("$BAM_OUT/$fullID.sorted.bam") 		|| print $logprint "<WARN>\t",timer(),"\tCan't delete $fullID\.sorted.bam: No such file!\n";
		unlink("$BAM_OUT/$fullID.sorted.bam.bai") 	|| print $logprint "<WARN>\t",timer(),"\tCan't delete $fullID\.sorted.bam\.bai: No such file!\n";
		# Renaming.
		move("$BAM_OUT/$fullID.nodup.bam","$BAM_OUT/$fullID.bam") 		|| die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
		move("$BAM_OUT/$fullID.nodup.bam.bai","$BAM_OUT/$fullID.bam.bai") 	|| die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
		# Finished.
		print $logprint "<INFO>\t",timer(),"\tFinished mapping for $fullID!\n";
	}
	undef(%input);
}


1;
