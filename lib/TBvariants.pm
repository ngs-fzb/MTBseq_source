#!/usr/bin/env perl

package TBvariants;

use strict;
use warnings;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION =  1.0.1;
@ISA     =  qw(Exporter);
@EXPORT  =  qw(tbvariants);

sub tbvariants {
   # get parameter and input from front-end.
   my $logprint            =  shift;
   my $VAR_dir             =  shift;
   my $POS_OUT             =  shift;
   my $CALL_OUT            =  shift;
   my $ref                 =  shift;
   my $refg                =  shift;
   my $micovf              =  shift;
   my $micovr              =  shift;
   my $miphred20           =  shift;
   my $mifreq              =  shift;
   my $resi_list_master    =  shift;
   my $int_regions         =  shift;
   my $all_vars            =  shift;
   my $snp_vars            =  shift;
   my $lowfreq_vars        =  shift;
   my @pos_files           =  @_;
   my $annotation          =  {};
   my $genes               =  {};
   my $variant_infos       =  {};
   my $resi_gene           =  {};
   # parse the genomic sequence for determining substitutions.
   print $logprint "<INFO>\t",timer(),"\tStart parsing $ref...\n";
   my $genome  =  parse_fasta($logprint,$VAR_dir,$ref);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
   # parse all necessary annotation information.
   print $logprint "<INFO>\t",timer(),"\tParsing $refg...\n";
   parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
   print $logprint "<INFO>\t",timer(),"\tStart parsing $resi_list_master and $int_regions...\n";
   parse_variant_infos($logprint,$variant_infos,$resi_gene,$resi_list_master,$int_regions);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $resi_list_master and $int_regions!\n";
   # start logic...
   foreach my $file (sort { $a cmp $b } @pos_files) {
      my $position_table   =  {};
      my $variants         =  {};
      my $statistics       =  {};
      # file names for output.
      $file                =~ /^(.*)_position_table.*\.tab$/;
      my $id               =  $1;
      my $output_mode      =  $all_vars . $snp_vars . $lowfreq_vars;
      my $param_string     =  "cf" . "$micovf" . "_cr" . "$micovr" . "_fr" . "$mifreq" . "_ph" . "$miphred20" . "_outmode" . "$output_mode";
      my $variants_file    =  "$id" . "_position_variants_" . $param_string . ".tab";
      my $uncovered_file   =  "$id" . "_position_uncovered_" . $param_string . ".tab";
      # parse position tables, call variants, make statistics, and calculate exchanges.
      print $logprint "<INFO>\t",timer(),"\tStart parsing $file...\n";
      parse_position_table($logprint,$POS_OUT,$file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
      print $logprint "<INFO>\t",timer(),"\tFinished parsing $file!\n";
      print $logprint "<INFO>\t",timer(),"\tStart calling variants from $file...\n";
      call_variants($logprint,$position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
      $position_table      =  {};
      $statistics          =  {};
      print $logprint "<INFO>\t",timer(),"\tFinished calling variants from $file!\n";
      print $logprint "<INFO>\t",timer(),"\tPrinting:\n";
      print $logprint "<INFO>\t",timer(),"\t$variants_file\n";
      print $logprint "<INFO>\t",timer(),"\t$uncovered_file\n";
      print_variants($CALL_OUT,$variants,$variant_infos,$resi_gene,$variants_file,$uncovered_file,$genes);
      $variants            =  {};
      print $logprint "<INFO>\t",timer(),"\tPrinting finished!\n";
   }
   @pos_files        =  ();
   $annotation       =  {};
   $genes            =  {};
   $variant_infos    =  {};
   $resi_gene        =  {};
}

1;
