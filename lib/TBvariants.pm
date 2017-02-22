#!/usr/bin/env perl

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

package TBvariants;

use strict;
use warnings;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###                                                                                                             ###
### Description: This package creates a variant lists from your input samples					###
###                                                                                                             ###
### Input:  .gatk_position_table.tab                   								###
### Output: .gatk_position_variants_cf*_cr*_fr*_ph*_out_mode*.tab						###
###	    .gatk_position_statistics_cf*_cr*_fr*_ph*_out_mode*.tab						###
###	    .gatk_position_uncovered_cf*_cr*_fr*_ph*_out_mode*.tab                                              ###
###                                                                                                             ###
###################################################################################################################

$VERSION	=	1.00;
@ISA		= 	qw(Exporter);
@EXPORT		= 	qw(tbvariants);

my $tbvariants;

sub tbvariants {
	# Switches autoflush for direct printing on.
	$|				=	1;
	# Get parameter and input from front-end.
	my $VAR_dir			=	shift;
	my $POS_OUT			=      	shift;
	my $CALL_OUT                 	=      	shift;
        my $ref                         =       shift;
        my $refg                        =       shift;	
	my $micovf			=	shift;
	my $micovr			=	shift;
	my $miphred20			=	shift;
	my $mifreq			=	shift;
	my $resi_list_master		=	shift;
	my $int_regions			=	shift;
	my $output_mode			=	shift;
	my @pos_files			=	@_;
    	my $annotation			=	{};
    	my $genes			=	{};
	my $variant_infos		=	{};
	my $resi_gene			=	{};
	# Parse the genomic sequence for determining substitutions.
	print  "<INFO>\t",timer(),"\tParsing $ref...\n";
	my $genome			=	parse_fasta($VAR_dir,$ref);
    	print  "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
	# Parse all necessary annotation information.
    	print  "<INFO>\t",timer(),"\tParsing $refg...\n";
	parse_annotation($VAR_dir,$genes,$annotation,$refg);
    	print  "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
	print  "<INFO>\t",timer(),"\tParsing $resi_list_master and $int_regions...\n";
	parse_variant_infos($variant_infos,$resi_gene,$resi_list_master,$int_regions);
	print  "<INFO>\t",timer(),"\tFinished parsing $resi_list_master and $int_regions!\n";
	# Start logic...
	foreach my $file(sort { $a cmp $b } @pos_files) {
		my $position_table		=	{};
		my $variants			=	{};
		my $statistics			=	{};
		print  "<INFO>\t",timer(),"\tStart variant calling from $file...\n";
		print  "<INFO>\t",timer(),"\tStart parsing $file...\n";
		# File names for output.
		$file				=~	/^(.*)_position_table.*\.tab$/;
		my $id 				=	$1;
		my $param_string		=	"cf" . "$micovf" . "_cr" . "$micovr" . "_fr" . "$mifreq" . "_ph" . "$miphred20" . "_outmode" . "$output_mode";
		my $statistics_file		=	"$id" . "_position_statistics_" . $param_string . ".tab";
		my $variants_file		=	"$id" . "_position_variants_" . $param_string . ".tab";
		my $uncovered_file		=	"$id" . "_position_uncovered_" . $param_string . ".tab";
		# Do the main work, meaning parse the position tables, call variants, make statistics, and calculate exchanges.
		parse_position_table($POS_OUT,$file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
		print  "<INFO>\t",timer(),"\tFinished parsing $file!\n";
		print  "<INFO>\t",timer(),"\tStart calling variants for $id...\n";
		call_variants($position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$output_mode);
		$position_table			=	{};
		print  "<INFO>\t",timer(),"\tFinished calling variants for $id!\n";
		print  "<INFO>\t",timer(),"\tPrinting:\n";
		print  "<INFO>\t",timer(),"\t$statistics_file\n";
		print  "<INFO>\t",timer(),"\t$variants_file\n";
		print  "<INFO>\t",timer(),"\t$uncovered_file\n";
		print_variants($CALL_OUT,$variants,$variant_infos,$resi_gene,$variants_file,$uncovered_file);
		$variants			=	{};
                print_stats($CALL_OUT,$statistics,$statistics_file,$id);
		$statistics			=	{};
		print  "<INFO>\t",timer(),"\tPrinting finished!\n";
		print  "<INFO>\t",timer(),"\tFinished variant calling from $file!\n";
	}
	$annotation				=	{};
	$genes					=	{};
	$variant_infos				=	{};
	$resi_gene				=	{};
}


1;
