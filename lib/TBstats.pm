#!/usr/bin/perl -w

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

package TBstats;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###                                                                                                             ###
### Description: This package use Samtools flagst for mapping statistics on the reference			###
###                                                                                                             ###
### Input:  .bam                                                 	                                      	###
### Output: mapping_statistics.txt                                                                              ###
###                                                                                                             ###
###################################################################################################################

$VERSION	=	1.00;
@ISA		=	qw(Exporter);
@EXPORT		=	qw(tbstats);


sub tbstats {
	# Switches autoflush on.
	$|			=	1;
	# Get parameter and input.
	my $W_dir		=	shift;
    	my $VAR_dir		=	shift;
	my $SAMTOOLS_dir        =       shift;
	my $BAM_OUT		=	shift;
	my $ref			=	shift;
	my $date_string         =       shift;
	my @bam_files		=	@_;
	my $stats_file		=	"mapping_statistics.tab";
 	my $input		=	{};
	my $statistics		=	{};
	my %check_up;
	my $genome		=	parse_fasta($VAR_dir,$ref);
	if(-f "$BAM_OUT/$stats_file") {
		open(IN,"$BAM_OUT/$stats_file");
		<IN>;
		while(<IN>) {
			next if($_ =~ /^\s/);
			my $line	=	$_;
			$line		=~      s/\015?\012?$//;
			my @fields	=	split(/\t/,$line);
			$check_up{$fields[1].$fields[2].$fields[3].$fields[4]} = 1;
		}
		close(IN);
	}
	# Start logic...
	foreach my $file(sort { $a cmp $b } @bam_files) {
		next unless (-f "$BAM_OUT/$file");
    		my @file_name		=	split(/_/,$file);
		my $sampleID		=	$file_name[0];
		my $libID		=	$file_name[1];
    		my $source		=	$file_name[2];
    		my $date		=	$file_name[3];
    		my $length		=	$file_name[4];
    		$length			=~	s/\.bam$//;
    		my $fullID		=	join("_",($sampleID,$libID,$source,$date,$length));
		if(exists $check_up{"\'$file_name[0]"."\'$file_name[1]"."\'$file_name[2]"."\'$file_name[3]"}) {
			print  "<INFO>\t",timer(),"\tSkipping, statistics calculation for $fullID. Is already existing!\n";
			next;
		}
		push(@{$input->{$fullID}},$file);
	}
	foreach my $fullID(sort { $a cmp $b } keys(%$input)) {
		my @bams		=	@{$input->{$fullID}};
    		if(scalar(@bams > 1)) {
			print  "<WARN>\t",timer(),"\tWarning, there is more than one .bam file for $fullID!\n";
			print  "<WARN>\t",timer(),"\tFiles will not be processed! check what went wrong with the naming scheme!\n";
			next;
		}
		print  "<INFO>\t",timer(),"\tStart using Samtools for statistics of BWA mapping for $fullID...\n";
		my @file_name		=	split("_",$bams[0]);
		my $sampleID		=	$file_name[0];
    		my $libID		=	$file_name[1];
    		my $source		=	$file_name[2];
    		my $date		=	$file_name[3];
    		my $length		=	$file_name[4];
    		$length			=~	s/bp\.bam$//;
		my $content		=	qx/$SAMTOOLS_dir\/samtools flagstat $BAM_OUT\/$fullID\.bam/;	
		die "$BAM_OUT/$fullID does not exist" if(!$content);
		my @lines		=	split(/\n/,$content);
		$lines[0]		=~	s/^(\d+)\s.*/$1/;
		$lines[4]		=~	s/^(\d+)\s.*/$1/;
		# Calculate theoretical coverage.
		my $tcov 		= 	($lines[0]*$length)/length($genome);
		$tcov			=	sprintf("%.2f",$tcov);
		# Calculate real coverage.
		my $rcov 		= 	($lines[4]*$length)/length($genome);
		$rcov    		= 	sprintf("%.2f",$rcov);
		# Calculate reads mapped in percent.
		my $relmap		=	0;
		$relmap 		=	($lines[4]/$lines[0])*100 if($lines[0] != 0);
		$relmap 		=  	sprintf("%.2f",$relmap);
		# Calculate a ratio between theoretical coverage and real coverage.
		my $ratio		=	0;
		$ratio			=	($tcov/$rcov) if($rcov != 0);
		$ratio   		=	sprintf("%.2f",$ratio);
		# Calculate the log2 of the ratio in order to be able to compare large differences.
		my $lratio		=	0;
		$lratio			=	log2(($tcov/$rcov)) if( $rcov != 0 );
		$lratio			=	sprintf("%.2f",$lratio);
		# Save everything into the hash.
		$statistics->{$source}->{$date}		.=	"\'$date_string"."\t"."\'$sampleID"."\t"."\'$libID"."\t"."\'$source"."\t"."\'$date"."\t"."\'$lines[0]"."\t"."\'$lines[4] \($relmap\)"."\t"."\'$tcov"."\t"."\'$rcov"."\t"."\'$ratio"."\t"."\'$lratio"."\n";
		print  "<INFO>\t",timer(),"\tFinished using Samtools for statistics of BWA mapping for $fullID!\n";
	}
	# Parse the hash.
    	print  "<INFO>\t",timer(),"\tStart printing statistics...\n";
	foreach my $source(sort { $a cmp $b } keys %$statistics) {
		foreach my $date(sort { $a cmp $b } keys %{$statistics->{$source}}) {
			unless(-f "$BAM_OUT/$stats_file") {
				open(OUT,">$BAM_OUT/$stats_file") || die "<INFO>\t",timer(),"\tCan't create $stats_file: $!\n";
				print OUT "Date\tSampleID\tLibraryID\tSource\tRun\tTotal_Reads\tMapped_Reads_(%)\tTheoretical_Coverage\tReal_Coverage\t(Theoretical/Real)\tlog2(Theoretical/Real)\n";	
				close(OUT);
			}
			open(OUT,">>$BAM_OUT/$stats_file") || die "<INFO>\t",timer(),"\tCan't create $stats_file: $!\n";
			print OUT $statistics->{$source}->{$date};
			close(OUT);
		}
	}
	# Finished.
   	print  "<INFO>\t",timer(),"\tFinished printing statistics!\n";
}


1;
