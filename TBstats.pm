#!/usr/bin/perl -w

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
### Output: Mapping_and_Variant_Statistics.txt                                                                  ###
###                                                                                                             ###
###################################################################################################################

$VERSION	=	1.10;
@ISA		=	qw(Exporter);
@EXPORT		=	qw(tbstats);


sub tbstats {
	# Get parameter and input.
	my $logprint		=	shift;
	my $W_dir		=	shift;
    	my $VAR_dir		=	shift;
	my $SAMTOOLS_dir        =       shift;
	my $BAM_OUT		=	shift;
	my $POS_OUT		=	shift;
	my $STATS_OUT		=	shift;
	my $ref			=	shift;
	my $refg		=	shift;
	my $micovf		=	shift;
	my $micovr		=	shift;
	my $miphred20		=	shift;
	my $mifreq		=	shift;
	my $all_vars		=	shift;
	my $snp_vars		=	shift;
	my $lowfreq_vars	=	shift;
	my $date_string         =       shift;
	my @bam_files		=	@_;
	my $stats_file		=	"Mapping_and_Variant_Statistics.tab";
	my $genes		=	{};
	my $annotation		=	{};
	my %check_up;
	my $genome		=	parse_fasta($logprint,$VAR_dir,$ref);
	# Parsing Annotation
	print $logprint "<INFO>\t",timer(),"\tParsing $refg...\n";
	parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
	print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
	if(-f "$STATS_OUT/$stats_file") {
		open(IN,"$STATS_OUT/$stats_file");
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
	foreach my $file_bam(sort { $a cmp $b } @bam_files) {
		next unless (-f "$BAM_OUT/$file_bam");
    		my @file_name		=	split(/_/,$file_bam);
		my $sampleID		=	$file_name[0];
		my $libID		=	$file_name[1];
    		my $source		=	$file_name[2];
    		my $date		=	$file_name[3];
    		my $length		=	$file_name[4];
    		$length			=~	s/(\d+).*$/$1/;
    		my $fullID		=	join("_",($sampleID,$libID,$source,$date,$length));
		if(exists $check_up{"\'$file_name[0]"."\'$file_name[1]"."\'$file_name[2]"."\'$file_name[3]"}) {
			print $logprint "<INFO>\t",timer(),"\tSkipping, statistics calculation for $fullID. Statisitcs already existing!\n";
			next;
		}
		print $logprint "<INFO>\t",timer(),"\tStart using Samtools for statistics of BWA mapping for $file_bam...\n";
		my $content		=	qx/$SAMTOOLS_dir\/samtools flagstat $BAM_OUT\/$file_bam/;	
		die print $logprint "$BAM_OUT/$file_bam does not exist" if(!$content);
		my @lines		=	split(/\n/,$content);
		$lines[0]		=~	s/^(\d+)\s.*/$1/;
		$lines[4]		=~	s/^(\d+)\s.*/$1/;
		# Calculate theoretical coverage.
		my $tcov 		= 	($lines[0] * $length) / length($genome);
		$tcov			=	sprintf("%.2f",$tcov);
		# Calculate real coverage.
		my $rcov 		= 	($lines[4] * $length) / length($genome);
		$rcov    		= 	sprintf("%.2f",$rcov);
		# Calculate reads mapped in percent.
		my $relmap		=	0;
		$relmap 		=	($lines[4] / $lines[0]) * 100 if($lines[0] != 0);
		$relmap 		=  	sprintf("%.2f",$relmap);
		# Calculate a ratio between theoretical coverage and real coverage.
		my $ratio		=	0;
		$ratio			=	($tcov / $rcov) if($rcov != 0);
		$ratio   		=	sprintf("%.2f",$ratio);
		# Calculate the log2 of the ratio in order to be able to compare large differences.
		my $lratio		=	0;
		$lratio			=	log2(($tcov / $rcov)) if($rcov != 0);
		$lratio			=	sprintf("%.2f",$lratio);
		print $logprint "<INFO>\t",timer(),"\tFinished using Samtools for statistics of BWA mapping for $file_bam!\n";
		my $pos_file		=	$fullID."\.gatk_position_table\.tab";
		print $logprint "<INFO>\t",timer(),"\tStart parsing position list of $pos_file...\n";
		# Get the variant statistics.
		my $position_table      =       {};
		my $variants            =       {};
		my $statistics          =       {};
		parse_position_table($logprint,$POS_OUT,$pos_file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
		print $logprint "<INFO>\t",timer(),"\tFinished parsing position list of $pos_file!\n";
		print $logprint "<INFO>\t",timer(),"\tStart fetching variant statistics from $pos_file...\n";
		call_variants($logprint,$position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
		$position_table		=	{};
		$variants		=	{};
		print $logprint "<INFO>\t",timer(),"\tFinished fetching variant statistics from $pos_file...\n";
		print $logprint "<INFO>\t",timer(),"\tStart preparing statistics for $fullID...\n";
		my $result              =       "'$date_string\t'$sampleID\t'$libID\t'$source\t'$date\t'$lines[0]\t'$lines[4] ($relmap)\t'$tcov\t'$rcov\t'$ratio\t'$lratio\t";
		$result = $result.prepare_stats($statistics);
		$statistics		=	{};
		print $logprint "<INFO>\t",timer(),"\tFinished preparing statistics for $fullID...\n";
		# Print statistics.
		print $logprint "<INFO>\t",timer(),"\tStart printing statistics...\n";
		unless(-f "$STATS_OUT/$stats_file") {
			open(OUT,">$STATS_OUT/$stats_file") || die print $logprint "<INFO>\t",timer(),"\tCan't create $stats_file: $!\n";
			my $header 	=	"Date\tSampleID\tLibraryID\tSource\tRun\tTotal Reads\tMapped Reads (%)\tTheoretical Coverage\tReal Coverage\t(Theoretical/Real)\tlog2(Theoretical/Real)\t";
			$header 	= 	$header."Genome Size\tGenome GC\t(Any) Total Bases (%)\t(Any) GC-Content\t(Any)Coverage mean/median\t(Unambiguous) Total Bases (%)\t(Unambiguous) GC-Content\t(Unambiguous) Coverage mean/median\t";
			$header		=	$header."SNPs\tDeletions\tInsertions\tUncovered\tSubstitutions (including Stop Codons)\n";
			print OUT $header;
			close(OUT);
		}
		open(OUT,">>$STATS_OUT/$stats_file") || die print $logprint "<INFO>\t",timer(),"\tCan't create $stats_file: $!\n";
		print OUT $result;
		close(OUT);
	}
	$genes		=	{};
	$annotation	=	{};
	# Finished.
   	print $logprint "<INFO>\t",timer(),"\tFinished printing statistics!\n";
}


1;
