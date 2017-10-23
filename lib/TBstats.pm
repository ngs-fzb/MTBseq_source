#!/usr/bin/perl -w

# tabstop is set to 8.

package TBstats;

use strict;
use warnings;
use diagnostics;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION	=	1.0.0;
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
	# Parsing Genome.
	print $logprint "<INFO>\t",timer(),"\tStart parsing $ref...\n";
	my $genome		=	parse_fasta($logprint,$VAR_dir,$ref);
	print $logprint "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
	# Parsing Annotation.
	print $logprint "<INFO>\t",timer(),"\tStart parsing $refg...\n";
	parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
	print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
	# Save existing statistics.
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
	@bam_files = sort { $a cmp $b } @bam_files;
	print "FILES ciao" . join("@",@bam_files)."\n";
	foreach my $file (@bam_files) {
		print "FILE: $file\n";
		
			next unless (-f "$BAM_OUT/$file");
    			my @file_name		=	split(/_/,$file);
			my $sampleID		=	shift(@file_name);
		print "SAMPLE: $sampleID\n";
			my $libID		=	shift(@file_name);
			$libID		=~	s/\.bam//;
    		print "LIBID: $libID\n";
			my $file_mod		=	join("_",@file_name);
		print "FILE_MOD: $file_mod\n";
			my $date		=	$1;
			my $source		=	$file_mod;
    			my $fullID		=	join("_",($sampleID,$libID));
			print "FULLID: $fullID\n";
			if($source ne ""){
			print "SOURCE: $source\n";
			my $source_new = substr($source,0,(length($source)-4));
			print "SOURCE_NEW: $source_new\n";
			$fullID.="_".$source_new;
			print "FULLID: $fullID\n";
		}
	
		# Check if sttistics for sample already exist.
		if(exists $check_up{"\'$sampleID"."\'$libID"}) {
			print $logprint "<INFO>\t",timer(),"\tSkipping, statistics calculation for $fullID. Statisitcs already existing!\n";
			next;
		}
		print $logprint "<INFO>\t",timer(),"\tStart using Samtools for BWA mapping statistics of $file...\n";
		my $content		=	qx/$SAMTOOLS_dir\/samtools flagstat $BAM_OUT\/$file/;	
		die print $logprint "$BAM_OUT/$file does not exist, TBstats line: ", __LINE__ , " \n" if(!$content);
		my @lines		=	split(/\n/,$content);
		$lines[0]		=~	s/^(\d+)\s.*/$1/;
		$lines[4]		=~	s/^(\d+)\s.*/$1/;
		# Calculate reads mapped in percent.
		my $relmap		=	0;
		$relmap 		=	($lines[4] / $lines[0]) * 100 if($lines[0] != 0);
		$relmap 		=  	sprintf("%.2f",$relmap);
		print $logprint "<INFO>\t",timer(),"\tFinished using Samtools for BWA mapping statistics of $file!\n";
		#my $pos_file_cor	=	substr($file,0,(length($file)-4));
		my $pos_file		=	$fullID."\.gatk_position_table\.tab";
		print $logprint "<INFO>\t",timer(),"\tStart parsing position list of $pos_file...\n";
		# Get variant statistics.
		my $position_table      =       {};
		my $variants            =       {};
		my $statistics          =       {};
		parse_position_table($logprint,$POS_OUT,$pos_file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
		print $logprint "<INFO>\t",timer(),"\tFinished parsing position list of $pos_file!\n";
		print $logprint "<INFO>\t",timer(),"\tStart fetching variant statistics from $pos_file...\n";
		call_variants($logprint,$position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
		$position_table		=	{};
		$variants		=	{};
		print $logprint "<INFO>\t",timer(),"\tFinished fetching variant statistics from $pos_file!\n";
		print $logprint "<INFO>\t",timer(),"\tStart preparing statistics for $fullID...\n";
		my $result              =       "'$date_string\t'$sampleID\t'$libID\t'$lines[0]\t'$lines[4]\t'$relmap\t";
		$result 		=	$result.prepare_stats($statistics);
		$statistics		=	{};
		print $logprint "<INFO>\t",timer(),"\tFinished preparing statistics for $fullID!\n";
		# Print statistics.
		print $logprint "<INFO>\t",timer(),"\tStart printing statistics into $stats_file...\n";
		unless(-f "$STATS_OUT/$stats_file") {
			open(OUT,">$STATS_OUT/$stats_file") || die print $logprint "<INFO>\t",timer(),"\tCan't create $stats_file: TBstats line: ", __LINE__ , " \n";
			my $header 	=	"Date\tSampleID\tLibraryID\tTotal Reads\tMapped Reads\t% Mapped Reads\t";
			$header 	= 	$header."Genome Size\tGenome GC\t(Any) Total Bases\t% (Any) Total Bases\t(Any) GC-Content\t(Any) Coverage mean\t(Any) Coverage median\t(Unambiguous) Total Bases\t% (Unambiguous) Total Bases\t(Unambiguous) GC-Content\t(Unambiguous) Coverage mean\t(Unambiguous) Coverage median\t";
			$header		=	$header."SNPs\tDeletions\tInsertions\tUncovered\tSubstitutions (Including Stop Codons)\n";
			print OUT $header;
			close(OUT);
		}
		open(OUT,">>$STATS_OUT/$stats_file") || die print $logprint "<INFO>\t",timer(),"\tCan't create $stats_file: TBstats line: ", __LINE__ , " \n";
		print OUT $result;
		close(OUT);
	}
	@bam_files      =       ();
	$genes		=	{};
	$annotation	=	{};
	undef(%check_up);
	# Finished.
   	print $logprint "<INFO>\t",timer(),"\tFinished printing statistics into $stats_file!\n";
}


1;
