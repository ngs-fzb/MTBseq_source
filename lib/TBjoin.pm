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

package TBjoin;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###                                                                                                             ###
### Description: This package creates a joint variant list from input samples					###
###                                                                                                             ###
### Input:  .gatk_position_table.tab, .gatk_position_variants_cf*_cr*_fr*_ph*_outmode*.tab                      ###
### Output: joint_cf*_cr*_fr*_ph*_outmode*_samples*.tab                                                         ###
###                                                                                                             ###
###################################################################################################################

$VERSION	=	1.10;
@ISA		=	qw(Exporter);
@EXPORT		=	qw(tbjoin);


sub tbjoin {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $VAR_dir		=	shift;
	my $POS_OUT		=	shift;
	my $CALL_OUT		=	shift;
	my $JOIN_OUT		=	shift;
	my $group_name		=	shift;
	my $ref                 =       shift;
        my $refg                =       shift;
	my $micovf		=      	shift;
	my $micovr		=      	shift;
	my $miphred20		=      	shift;
	my $mifreq		=      	shift;
	my $all_vars		=	shift;
	$all_vars		=	1;
	my $snp_vars		=	shift;
	my $lowfreq_vars	=	shift;
	my @var_files		=	@_;
	my $annotation		=	{};
	my $genes		=	{};
	my $var_positions	=	{};
	my $strain		=	{};
	my $position_stats	=	{};
	my $param_string	=	"_cf" . "$micovf" . "_cr" . "$micovr" . "_fr" . "$mifreq" . "_ph" . "$miphred20";
	my $join_file		=	"$group_name" . "_joint" . "$param_string" . "_samples" . scalar(@var_files) . ".tab";
	my $breadth_file	=	"$group_name" . "_joint" . "$param_string" . "_samples" . scalar(@var_files) . ".log";
	my @ids;
	# Get the genomic sequence for determining substitutions.
	print $logprint "<INFO>\t",timer(),"\tParsing $ref...\n";
	my $genome		=	parse_fasta($logprint,$VAR_dir,$ref);
    	print $logprint "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
	# Prepare info for coverage_breadth output.
	for(my $i = 1; $i <= length($genome); $i++) {
		$position_stats->{$i}->{0}->{unambigous}	=	0;
		$position_stats->{$i}->{0}->{any}        	=       0;
		$position_stats->{$i}->{0}->{nothing}        	=       0;
	}
	# Get all necessary annotation information.
    	print $logprint "<INFO>\t",timer(),"\tParsing $refg...\n";	
	parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
    	print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg...\n";
	# Parse all variant files to $positions_hash and @ids
	print $logprint "<INFO>\t",timer(),"\tParsing variant files...\n";
	foreach my $file(sort { $a cmp $b } @var_files) {
		$file			=~	/^(.*)\.gatk_position_variants_.*\.tab$/;
		my $id			=	$1;
		parse_variants($logprint,$CALL_OUT,$file,$id,$var_positions,$strain,$micovf,$micovr,$mifreq,$miphred20);
		push(@ids, $id);
	}
	print $logprint "<INFO>\t",timer(),"\tFinished parsing variant files!\n";
    	print $logprint "<INFO>\t",timer(),"\tPrinting joint variant file scaffold...\n";
	print_joint_table_scaffold($logprint,$JOIN_OUT,$join_file,$var_positions,$annotation,$genes,@ids);
    	print $logprint "<INFO>\t",timer(),"\tFinished printing joint variant file scaffold!\n";
	# Start logic..
	print $logprint "<INFO>\t",timer(),"\tParsing position lists, extend called variants and complete joint variant list...\n";
	foreach my $id (@ids) {
		my $position_file		=	$id.".gatk_position_table.tab";
		next unless(-f "$POS_OUT/$position_file");
		# Parse position_table.
		my $position_table		=	{};
		parse_position_table($logprint,$POS_OUT,$position_file,$micovf,$micovr,$miphred20,$mifreq,$position_table,$position_stats);
		# Skip every position not included in $var_positions.
		my $joint_position_table	=	{};
		foreach my $pos (keys %$var_positions) {
			foreach my $insertion_index (keys %{$var_positions->{$pos}}) {
				if($insertion_index != 0) {
					$joint_position_table->{$pos}->{$insertion_index}	=	$position_table->{$pos}->{$insertion_index} if(exists $position_table->{$pos}->{$insertion_index});
					$joint_position_table->{$pos}->{"0"}       		=       $position_table->{$pos}->{"0"};
				}
				else {
					$joint_position_table->{$pos}->{$insertion_index}       =       $position_table->{$pos}->{$insertion_index};
				}
			}
		}
		$position_table			=	{};
		my $variants			=	{};
		my $statistics			=	{};
		call_variants($logprint,$joint_position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
		$joint_position_table		=	{};
		$statistics			=	{};
		# This already prints the information for the dataset to save RAM.
		print_joint_table($logprint,$JOIN_OUT,$join_file,$id,$strain,$variants);
	}
	$annotation		=	{};
        $genes			=	{};
        $var_positions		=	{};
        $strain			=	{};
	print $logprint "<INFO>\t",timer(),"\tStart printing coverage breadth output...\n";
	print_position_stats($JOIN_OUT,$breadth_file,$genome,$position_stats,@ids);
	print $logprint "<INFO>\t",timer(),"\tFinished printing coverage breadth output!\n";
	$position_stats		=	{};
}


1;
