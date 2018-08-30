#!/usr/bin/perl

package TBrefine2;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION =  1.0.0;
@ISA     =  qw(Exporter);
@EXPORT  =  qw(tbrefine2);

sub tbrefine2 {
   # get parameter and input from front-end.
   my $logprint      =  shift;
   my $W_dir         =  shift;
   my $VAR_dir       =  shift;
   my $PICARD_dir    =  shift;
   my $GATK_dir      =  shift;
   my $BAM_OUT       =  shift;
   my $GATK_OUT      =  shift;
   my $ref           =  shift;
   my $basecalib     =  shift;
   my $threads       =  shift;
   my @bam_files     =  @_;
   my %input;
   # start logic...
   foreach my $file (sort {$a cmp $b } @bam_files) {
      my @file_name     =  split(/_/,$file);
      my $sampleID      =  shift(@file_name);
      my $libID         =  shift(@file_name);
      $libID            =~ s/\.bam//;
      my $file_mod      =  join("_",@file_name);
      my $source        =  $file_mod;
      my $fullID        =  join("_",($sampleID,$libID));
      if($source ne "") {
         my $source_new =  substr($source,0,(length($source)-4));
         $fullID.="_".$source_new;
      }
      push(@{$input{$fullID}},$file);
   }
   @bam_files           =  ();
   foreach my $fullID (sort { $a cmp $b } keys %input) {
      my @bams          =  @{$input{$fullID}};
      if(scalar(@bams > 1)) {
         print $logprint "<WARN>\t",timer(),"\tSkipping $fullID, more than one .bam file for $fullID!\n";
         next;
      }
      my @file_name     =  split(/_/,$fullID);
      my $sampleID      =  shift(@file_name);
      my $libID         =  shift(@file_name);
      my $file_mod      =  join("_",@file_name);
      my $source        =  $file_mod;
      my $date          =  $1;
      print $logprint "<INFO>\t",timer(),"\tUpdating log file for $fullID...\n";
      my $old_logfile   =  $fullID.".bamlog";
      my $logfile       =  $fullID.".gatk.bamlog";
      unlink("$GATK_OUT/$logfile");
      if(-f "$BAM_OUT/$old_logfile") {
         cat($logprint,"$BAM_OUT/$old_logfile","$GATK_OUT/$logfile") || die print $logprint "<ERROR>\t",timer(),"\tcat failed: TBrefine.pm line: ", __LINE__ , " \n";
      }
      my $dict          =  $ref;
      $dict             =~ s/\.fasta/.dict/;
      # create a dictionary with picard tools, if it isn't created already.
      unless(-f "$VAR_dir/$dict") {
         print $logprint "<INFO>\t",timer(),"\tStart using Picard Tools for creating a dictionary of $ref...\n";
         print $logprint "<INFO>\t",timer(),"\tjava -jar $PICARD_dir/picard.jar CreateSequenceDictionary R=$VAR_dir/$ref O=$VAR_dir/$dict 2>> $GATK_OUT/$logfile\n";
         system("java -jar $PICARD_dir/picard.jar CreateSequenceDictionary R=$VAR_dir/$ref O=$VAR_dir/$dict 2>> $GATK_OUT/$logfile");
         print $logprint "<INFO>\t",timer(),"\tFinished using Picard Tools for creating a dictionary of $ref!\n";
      }
      ## use RealignerTargetCreator from GATK. OUTDATED 13-07-2018
      #print $logprint "<INFO>\t",timer(),"\tStart using GATK RealignerTargetCreator for $fullID...\n";
      #print $logprint "<INFO>\t",timer(),"\tjava -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type RealignerTargetCreator --reference_sequence $VAR_dir/$ref --input_file $BAM_OUT/$fullID.bam --downsample_to_coverage 10000 --num_threads $threads --out $GATK_OUT/$fullID.gatk.intervals 2>> $GATK_OUT/$logfile\n";
      #system("java -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type RealignerTargetCreator --reference_sequence $VAR_dir/$ref --input_file $BAM_OUT/$fullID.bam --downsample_to_coverage 10000 --num_threads $threads --out $GATK_OUT/$fullID.gatk.intervals 2>> $GATK_OUT/$logfile");
      #print $logprint "<INFO>\t",timer(),"\tFinished using GATK RealignerTargetCreator for $fullID!\n";
      ## use IndelRealigner from GATK.
      #print $logprint "<INFO>\t",timer(),"\tStart using GATK IndelRealigner for $fullID...\n";
     # print $logprint "<INFO>\t",timer(),"\tjava -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type IndelRealigner --reference_sequence $VAR_dir/$ref --input_file $BAM_OUT/$fullID.bam --defaultBaseQualities 12 --targetIntervals $GATK_OUT/$fullID.gatk.intervals --noOriginalAlignmentTags --out $GATK_OUT/$fullID.realigned.bam 2>> $GATK_OUT/$logfile\n";
      #system("java -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type IndelRealigner --reference_sequence $VAR_dir/$ref --input_file $BAM_OUT/$fullID.bam --defaultBaseQualities 12 --targetIntervals $GATK_OUT/$fullID.gatk.intervals --noOriginalAlignmentTags --out $GATK_OUT/$fullID.realigned.bam 2>> $GATK_OUT/$logfile");
      #print $logprint "<INFO>\t",timer(),"\tFinished using GATK IndelRealigner for $fullID!\n";
      ## if $basecalib is not defined, skip the next parts.
      if($basecalib eq 'NONE') {
         print $logprint "<INFO>\t",timer(),"\tSkipping GATK BaseRecalibrator! No calibration file specified!\n";
         move("$BAM_OUT/$fullID.bam","$GATK_OUT/$fullID.gatk.bam") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: TBrefine.pm line: ", __LINE__ , " \n";
         move("$BAM_OUT/$fullID.bai","$GATK_OUT/$fullID.gatk.bai") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: TBrefine.pm line: ", __LINE__ , " \n";
         next;
      }
	  # use BaseRecalibrator from GATK.
      print $logprint "<INFO>\t",timer(),"\tStart using GATK4 BaseRecalibrator for $fullID...\n";
      print $logprint "<INFO>\t",timer(),"\t$GATK_dir/gatk BaseRecalibrator --reference $VAR_dir/$ref --input $BAM_OUT/$fullID.bam --known-sites $basecalib --maximum-cycle-value 600 --output $GATK_OUT/$fullID.gatk.recal.table 2>> $GATK_OUT/$logfile\n";
      system("$GATK_dir/gatk BaseRecalibrator --reference $VAR_dir/$ref --input $BAM_OUT/$fullID.bam --known-sites $basecalib --maximum-cycle-value 600 --output $GATK_OUT/$fullID.gatk.recal.table 2>> $GATK_OUT/$logfile");
      print $logprint "<INFO>\t",timer(),"\tFinished using GATK BaseRecalibrator for $fullID!\n";
      # use ApplyBQSR with GATK.
      print $logprint "<INFO>\t",timer(),"\tStart using GATK ApplyBQSR for $fullID...\n";
      print $logprint "<INFO>\t",timer(),"\t$GATK_dir/gatk ApplyBQSR --reference $VAR_dir/$ref --input $BAM_OUT/$fullID.bam -bqsr $GATK_OUT/$fullID.gatk.recal.table --output $GATK_OUT/$fullID.gatk.bam  2>> $GATK_OUT/$logfile\n";
      system("$GATK_dir/gatk ApplyBQSR --reference $VAR_dir/$ref --input $BAM_OUT/$fullID.bam -bqsr $GATK_OUT/$fullID.gatk.recal.table --output $GATK_OUT/$fullID.gatk.bam  2>> $GATK_OUT/$logfile");
      print $logprint "<INFO>\t",timer(),"\tFinished using GATK ApplyBQSR for $fullID!\n";
      # removing temporary files.
      #print $logprint "<INFO>\t",timer(),"\tRemoving temporary files...\n";
      #unlink("$GATK_OUT/$fullID.realigned.bam");
      #unlink("$GATK_OUT/$fullID.realigned.bai");
      # finished.
      print $logprint "<INFO>\t",timer(),"\tGATK refinement finished for $fullID!\n";
   }
}

1;
