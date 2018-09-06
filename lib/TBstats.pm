#!/usr/bin/env perl -w

package TBstats;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION =  1.0.0;
@ISA     =  qw(Exporter);
@EXPORT  =  qw(tbstats);

sub tbstats {
   # get parameter and input.
   my $logprint      =  shift;
   my $W_dir         =  shift;
   my $VAR_dir       =  shift;
   my $SAMTOOLS_dir  =  shift;
   my $SAMTOOLS_call =  shift;
   my $BAM_OUT       =  shift;
   my $POS_OUT       =  shift;
   my $STATS_OUT     =  shift;
   my $ref           =  shift;
   my $refg          =  shift;
   my $micovf        =  shift;
   my $micovr        =  shift;
   my $miphred20     =  shift;
   my $mifreq        =  shift;
   my $all_vars      =  shift;
   my $snp_vars      =  shift;
   my $lowfreq_vars  =  shift;
   my $date_string   =  shift;
   my @bam_files     =  @_;
   my $stats_file    =  "Mapping_and_Variant_Statistics.tab";
   my $genes         =  {};
   my $annotation    =  {};
   my %check_up;
   # parsing genome.
   print $logprint "<INFO>\t",timer(),"\tStart parsing $ref...\n";
   my $genome        =  parse_fasta($logprint,$VAR_dir,$ref);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
   # parsing annotation.
   print $logprint "<INFO>\t",timer(),"\tStart parsing $refg...\n";
   parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
   # save existing statistics.
   if(-f "$STATS_OUT/$stats_file") {
      open(IN,"$STATS_OUT/$stats_file");
      <IN>;
      while(<IN>) {
         next if($_ =~ /^\s/);
         my $line    =  $_;
         $line       =~ s/\015?\012?$//;
         my @fields  =  split(/\t/,$line);
         $check_up{$fields[1].$fields[2]} = 1;
      }
      close(IN);
   }
   # start logic...
   foreach my $file (sort { $a cmp $b } @bam_files) {
      next unless (-f "$BAM_OUT/$file");
      my @file_name     =  split(/_/,$file);
      my $sampleID      =  shift(@file_name);
      my $libID         =  shift(@file_name);
      $libID            =~ s/\.bam//;
      my $file_mod      =  join("_",@file_name);
      my $source        =  $file_mod;
      my $fullID        =  join("_",($sampleID,$libID));
      if($source ne "") {
         my $source_new =  substr($source,0,(length($source)-4));
         $fullID        .= "_" . $source_new;
      }
		# check if statistics for sample already exist.
      if(exists $check_up{"\'$sampleID"."\'$libID"}) {
         print $logprint "<INFO>\t",timer(),"\tSkipping, statistics calculation for $fullID Statisitcs already existing!\n";
         next;
      }
      print $logprint "<INFO>\t",timer(),"\tStart using Samtools for BWA mapping statistics of $file...\n";
      my $content    =  qx/$SAMTOOLS_call flagstat $BAM_OUT\/$file/;
      die print $logprint "$BAM_OUT/$file does not exist, TBstats.pm line: ", __LINE__ , " \n" if(!$content);
      my @lines      =  split(/\n/,$content);
      $lines[0]      =~ s/^(\d+)\s.*/$1/;
      $lines[4]      =~ s/^(\d+)\s.*/$1/;
      # calculate reads mapped in percent.
      my $relmap     =  0;
      $relmap        =  ($lines[4] / $lines[0]) * 100 if($lines[0] != 0);
      $relmap        =  sprintf("%.2f",$relmap);
      print $logprint "<INFO>\t",timer(),"\tFinished using Samtools for BWA mapping statistics of $file!\n";
      my $pos_file   =  $fullID . "\.gatk_position_table\.tab";
      print $logprint "<INFO>\t",timer(),"\tStart parsing position list of $pos_file...\n";
      # get variant statistics.
      my $position_table   =  {};
      my $variants         =  {};
      my $statistics       =  {};
      parse_position_table($logprint,$POS_OUT,$pos_file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
      print $logprint "<INFO>\t",timer(),"\tFinished parsing position list of $pos_file!\n";
      print $logprint "<INFO>\t",timer(),"\tStart fetching variant statistics from $pos_file...\n";
      call_variants($logprint,$position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
      $position_table      =  {};
      $variants            =  {};
      print $logprint "<INFO>\t",timer(),"\tFinished fetching variant statistics from $pos_file!\n";
      print $logprint "<INFO>\t",timer(),"\tStart preparing statistics for $fullID...\n";
      my $result           =  "'$date_string\t'$sampleID\t'$libID\t'$fullID\t'$lines[0]\t'$lines[4]\t'$relmap\t";
      $result              =  $result.prepare_stats($statistics);
      $statistics          =  {};
      print $logprint "<INFO>\t",timer(),"\tFinished preparing statistics for $fullID!\n";
      # print statistics.
      print $logprint "<INFO>\t",timer(),"\tStart printing statistics into $stats_file...\n";
      unless(-f "$STATS_OUT/$stats_file") {
         open(OUT,">$STATS_OUT/$stats_file") || die print $logprint "<INFO>\t",timer(),"\tCan't create $stats_file: TBstats.pm line: ", __LINE__ , " \n";
         my $header        =  "Date\tSampleID\tLibraryID\tFullID\tTotal Reads\tMapped Reads\t% Mapped Reads\t";
         $header           =  $header."Genome Size\tGenome GC\t(Any) Total Bases\t% (Any) Total Bases\t(Any) GC-Content\t(Any) Coverage mean\t(Any) Coverage median\t(Unambiguous) Total Bases\t% (Unambiguous) Total Bases\t(Unambiguous) GC-Content\t(Unambiguous) Coverage mean\t(Unambiguous) Coverage median\t";
         $header           =  $header."SNPs\tDeletions\tInsertions\tUncovered\tSubstitutions (Including Stop Codons)\n";
         print OUT $header;
         close(OUT);
      }
      open(OUT,">>$STATS_OUT/$stats_file") || die print $logprint "<INFO>\t",timer(),"\tCan't create $stats_file: TBstats line: ", __LINE__ , " \n";
      print OUT $result;
      close(OUT);
   }
   @bam_files     =  ();
   $genes         =  {};
   $annotation    =  {};
   undef(%check_up);
   # finished.
   print $logprint "<INFO>\t",timer(),"\tFinished printing statistics into $stats_file!\n";
}

1;
