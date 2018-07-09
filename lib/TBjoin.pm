#!/usr/bin/env perl

package TBjoin;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION    =  1.0.0;
@ISA        =  qw(Exporter);
@EXPORT     =  qw(tbjoin);

sub tbjoin {
   # get parameter and input from front-end.
   my $logprint         =  shift;
   my $VAR_dir          =  shift;
   my $POS_OUT          =  shift;
   my $CALL_OUT         =  shift;
   my $JOIN_OUT         =  shift;
   my $group_name       =  shift;
   my $ref              =  shift;
   my $refg             =  shift;
   my $micovf           =  shift;
   my $micovr           =  shift;
   my $miphred20        =  shift;
   my $mifreq           =  shift;
   my $all_vars         =  shift;
   $all_vars            =  1;
   my $snp_vars         =  shift;
   my $lowfreq_vars     =  shift;
   my @var_files        =  @_;
   my $annotation       =  {};
   my $genes            =  {};
   my $var_positions    =  {};
   my $strain           =  {};
   my $position_stats   =  {};
   my $param_string     =  "_cf" . "$micovf" . "_cr" . "$micovr" . "_fr" . "$mifreq" . "_ph" . "$miphred20";
   my $join_file        =  "$group_name" . "_joint" . "$param_string" . "_samples" . scalar(@var_files) . ".tab";
   my $breadth_file     =  "$group_name" . "_joint" . "$param_string" . "_samples" . scalar(@var_files) . ".log";
   my @ids;
   # get genomic sequence for determining substitutions.
   print $logprint "<INFO>\t",timer(),"\tStart parsing $ref...\n";
   my $genome  =  parse_fasta($logprint,$VAR_dir,$ref);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
   # prepare info for coverage_breadth output.
   for(my $i = 1; $i <= length($genome); $i++) {
      $position_stats->{$i}->{0}->{unambigous}  =  0;
      $position_stats->{$i}->{0}->{any}         =  0;
      $position_stats->{$i}->{0}->{nothing}     =  0;
   }
   # get all necessary annotation information.
   print $logprint "<INFO>\t",timer(),"\tStart parsing $refg...\n";
   parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg...\n";
   # parse all variant files to $positions_hash and @ids
   print $logprint "<INFO>\t",timer(),"\tParsing variant files...\n";
   foreach my $file(sort { $a cmp $b } @var_files) {
      $file    =~ /^(.*)\.gatk_position_variants_.*\.tab$/;
      my $id   =  $1;
      parse_variants($logprint,$CALL_OUT,$file,$id,$var_positions,$strain,$micovf,$micovr,$mifreq,$miphred20);
      push(@ids, $id);
   }
   print $logprint "<INFO>\t",timer(),"\tFinished parsing variant files!\n";
   print $logprint "<INFO>\t",timer(),"\tSart printing joint variant file scaffold...\n";
   print_joint_table_scaffold($logprint,$JOIN_OUT,$join_file,$var_positions,$annotation,$genes,@ids);
   print $logprint "<INFO>\t",timer(),"\tFinished printing joint variant file scaffold!\n";
   # start logic..
   print $logprint "<INFO>\t",timer(),"\tStart parsing position lists, extend called variants and complete joint variant list...\n";
   foreach my $id (@ids) {
      my $position_file          =  $id.".gatk_position_table.tab";
      next unless(-f "$POS_OUT/$position_file");
      # parse position_table.
      my $position_table         =  {};
      parse_position_table($logprint,$POS_OUT,$position_file,$micovf,$micovr,$miphred20,$mifreq,$position_table,$position_stats);
      # skip every position not included in $var_positions.
      my $joint_position_table   =  {};
      foreach my $pos (keys %$var_positions) {
         foreach my $insertion_index (keys %{$var_positions->{$pos}}) {
            if($insertion_index != 0) {
               $joint_position_table->{$pos}->{$insertion_index}  =  $position_table->{$pos}->{$insertion_index}  if(exists $position_table->{$pos}->{$insertion_index});
               $joint_position_table->{$pos}->{"0"}               =  $position_table->{$pos}->{"0"};
            }
            else {
               $joint_position_table->{$pos}->{$insertion_index}  =  $position_table->{$pos}->{$insertion_index};
            }
         }
      }
      $position_table         =  {};
      my $variants            =  {};
      my $statistics          =  {};
      call_variants($logprint,$joint_position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
      $joint_position_table   =  {};
      $statistics             =  {};
      # this print information for the dataset to save RAM.
      print_joint_table($logprint,$JOIN_OUT,$join_file,$id,$strain,$variants);
   }
   print $logprint "<INFO>\t",timer(),"\tFinished parsing position lists, extend called variants and complete joint variant list!\n";
   $annotation       =  {};
   $genes            =  {};
   $var_positions    =  {};
   $strain           =  {};
   print $logprint "<INFO>\t",timer(),"\tStart printing coverage breadth output...\n";
   print_position_stats($JOIN_OUT,$breadth_file,$genome,$position_stats,@ids);
   print $logprint "<INFO>\t",timer(),"\tFinished printing coverage breadth output!\n";
   $position_stats   =  {};
}

1;
