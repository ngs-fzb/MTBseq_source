#!/usr/bin/perl

package TBgroups;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION    =  1.0.0;
@ISA        =  qw(Exporter);
@EXPORT     =  qw(tbgroups);

sub tbgroups {
   # get parameter and input from front-end.
   my $logprint         =  shift;
   my $AMEND_OUT        =  shift;
   my $GROUPS_OUT       =  shift;
   my $distance         =  shift;
   my @joint_files      =  @_;
   my $strip_ids        =  1;
   # start logic...
   foreach my $joint_file (sort { $a cmp $b } @joint_files) {
      my $matrix_file      =  "";
      my $group_file       =  "";
      my $pivot_hash       =  {};
      my $position_info    =  {};
      my $distance_matrix  =  {};
      my $strains          =  [];
      $matrix_file         =  strip($joint_file,".matrix");
      $group_file          =  strip($joint_file,("_d" . $distance . ".groups"));
      print $logprint "<INFO>\t",timer(),"\tStart parsing $joint_file...\n";
      parse_amend_table($logprint,$AMEND_OUT,$joint_file,$pivot_hash,$strains,$position_info);
      print $logprint "<INFO>\t",timer(),"\tFinished parsing $joint_file!\n";
      print $logprint "<INFO>\t",timer(),"\tStart building distance matrix for $joint_file...\n";
      build_matrix($GROUPS_OUT,$matrix_file,$pivot_hash,$distance_matrix,$strains);
      print $logprint "<INFO>\t",timer(),"\tFinished building distance matrix for $joint_file!\n";
      print $logprint "<INFO>\t",timer(),"\tStart calling groups for $joint_file...\n";
      call_groups($GROUPS_OUT,$group_file,$distance,$strip_ids,$strains,$distance_matrix);
      print $logprint "<INFO>\t",timer(),"\tFinished calling groups for $joint_file!\n";
   }
}

1;
