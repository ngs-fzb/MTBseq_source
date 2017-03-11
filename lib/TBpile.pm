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

package TBpile;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###                                                                                                             ###
### Description: This package use GATK for variant calling. Currently, using the  HaplotypeCaller is under 	###
### construction. An .mpileup file is created in order to prepare an input for a precise home brewed variant 	###
### calling pipeline.                                        							###
###                                                                                                             ###
### Input:  .gatk.bam                                                                                           ###
### Output: [.gatk.vcf], [.gatk.vcflog], gatk.mpileup, .gatk.mpileuplog                                         ###
###														###
### [] = not produced at the moment.										###
###                                                                                                             ###
###################################################################################################################

$VERSION	=	1.11;
@ISA		=	qw(Exporter);
@EXPORT		=	qw(tbpile);


sub tbpile {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $VAR_dir		=	shift;
	my $SAMTOOLS_dir        =       shift;
	my $GATK_dir            =       shift;
	my $GATK_OUT		=	shift;
	my $MPILE_OUT		=	shift;
	my $ref                 =       shift;
	my $threads		=	shift;
	my @gbam_files		=	@_;
	# Start logic...
	foreach my $file (sort { $a cmp $b } @gbam_files) {
		my @file_name		=	split(/_/,$file);
		my $sampleID		= 	shift(@file_name);
		my $libID 		=	shift(@file_name);
		my $source		= 	shift(@file_name);
		my $date		= 	shift(@file_name);
		my $seqlength		=	shift(@file_name);
		$seqlength		=~	s/\.gatk\.bam//;
		my $fullID		=	join("_",($sampleID,$libID,$source,$date,$seqlength));
		print $logprint "<INFO>\t",timer(),"\tUpdating logfile for $fullID...\n";		
		my $basename		=	$file; 
		$basename		=~	s/\.bam$//;
		my $old_logfile		=	$basename.".bamlog";
		my $vcf_logfile		=	$basename.".vcflog";
		my $mpile_logfile	=	$basename.".mpileuplog";
		my $mpile_file		=	$basename.".mpileup";
		unlink("$GATK_OUT/$vcf_logfile")	|| print $logprint "<WARN>\t",timer(),"\tCan't delete $vcf_logfile: $!\n";
		unlink("$GATK_OUT/$mpile_logfile")	|| print $logprint "<WARN>\t",timer(),"\tCan't delete $mpile_logfile: $!\n";
		if(-f "$GATK_OUT/$old_logfile") {
			#cat($logprint,"$GATK_OUT/$old_logfile","$GATK_OUT/$vcf_logfile")	|| die "<ERROR>\t",timer(),"\tcat failed: $?\n";
			cat($logprint,"$GATK_OUT/$old_logfile","$GATK_OUT/$mpile_logfile") 	|| die "<ERROR>\t",timer(),"\tcat failed: $!\n";
		}

### UNDER CONSTRUCTION - Haplotype Caller and VCF Files ###
=head2
		
		# Variant calling
		print $logprint "<INFO>\t",timer(),"\tStart using GATK HaplotypeCaller for creating a raw callset of variants in $fullID...\n";
		system("java -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type HaplotypeCaller --reference_sequence $VAR_dir/$ref --input_file $GATK_OUT/$file --genotyping_mode DISCOVERY --standard_min_confidence_threshold_for_emitting 10 --standard_min_confidence_threshold_for_calling 30 --sample_ploidy 1 --num_cpu_threads_per_data_thread $threads --emitRefConfidence BP_RESOLUTION --output_mode EMIT_ALL_SITES --out $GATK_OUT/$fullID.g.vcf 2>> $GATK_OUT/$vcf_logfile");
		print $logprint "<INFO>\t",timer(),"\tFinished using GATK HaplotypeCaller for creating a raw callset of variants in $fullID!\n";
		# Genotyping
        	print $logprint "<INFO>\t",timer(),"\tStart using GATK GenotypeGVCFs for genotyping of $fullID...\n";
		system("java -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type GenotypeGVCFs --reference_sequence $VAR_dir/$ref --variant $GATK_OUT/$fullID.g.vcf --standard_min_confidence_threshold_for_emitting 10 --standard_min_confidence_threshold_for_calling 30 --sample_ploidy 1 --num_threads $threads --includeNonVariantSites --out $GATK_OUT/$fullID.gatk.unfiltered.vcf 2>> $GATK_OUT/$vcf_logfile");
        	print $logprint "<INFO>\t",timer(),"\tFinished using GATK GenotypeGVCFs for genotyping of $fullID!\n";
		# Filtering SNPs
        	print $logprint "<INFO>\t",timer(),"\tStart using GATK VariantFiltration with default filter recommendation for SNPs in $fullID...\n";
       		system("java -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type VariantFiltration --reference_sequence $VAR_dir/$ref --variant $GATK_OUT/$fullID.gatk.unfiltered.vcf --filterExpression \"QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"SNP_filter\" --out $GATK_OUT/$fullID.gatk.filtered.vcf 2>> $GATK_OUT/$vcf_logfile");
        	print $logprint "<INFO>\t",timer(),"\tFinished using GATK VariantFiltration with default filter recommendation for SNPs in $fullID!\n";
        	print $logprint "<INFO>\t",timer(),"\tStart using GATK VariantFiltration with default filter recommendation for InDels in $fullID...\n";
		system("java -jar $GATK_dir/GenomeAnalysisTK.jar --analysis_type VariantFiltration --reference_sequence $VAR_dir/$ref --variant $GATK_OUT/$fullID.gatk.filtered.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0\" --filterName \"InDel_filter\" --out $GATK_OUT/$fullID.gatk.filtered.vcf 2>> $GATK_OUT/$vcf_logfile");
        	print $logprint "<INFO>\t",timer(),"\tFinished using GATK VariantFiltration with default filter recommendation for InDels in $fullID!\n";

=cut	
    		
		print $logprint "<INFO>\t",timer(),"\tRemoving temporary files for $fullID...\n";
    		#unlink("$GATK_OUT/$fullID.g.vcf")	|| print $logprint "<WARN>\t",timer(),"\tCan't delete $fullID.g.vcf: No such file!\n";
   		#unlink("$GATK_OUT/$fullID.g.vcf.idx") 	|| print $logprint "<WARN>\t",timer(),"\tCan't delete $fullID.g.vcf.idx: No such file!\n";
        	# Create .mpileup files.
        	print $logprint "<INFO>\t",timer(),"\tStart using samtools for creating a .mpileup file for $fullID...\n";
        	print $logprint "<INFO>\t",timer(),"\t$SAMTOOLS_dir/samtools mpileup -B -A -f $VAR_dir/$ref $GATK_OUT/$file > $MPILE_OUT/$mpile_file 2>> $MPILE_OUT/$mpile_logfile\n";
		system("$SAMTOOLS_dir/samtools mpileup -B -A -f $VAR_dir/$ref $GATK_OUT/$file > $MPILE_OUT/$mpile_file 2>> $MPILE_OUT/$mpile_logfile");
        	print $logprint "<INFO>\t",timer(),"\tFinished using samtools for creating a .mpileup file for $fullID!\n";
		# Finished.
        	print $logprint "<INFO>\t",timer(),"\tGATK variant calling finished for $fullID!\n";
	}
	@gbam_files = ();
}


1;
