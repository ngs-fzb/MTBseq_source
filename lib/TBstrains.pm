#!/usr/bin/perl

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

package TBstrains;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###                                                                                                             ###
### Description: This package uses Position Tables and Variant Calling for determining the underling Lineage.	###
###                                                          							###
###                                                                                                             ###
### Input:      .position_table.tab                                                                             ###
### Output:     strain_classification.tab                                                                       ###
###                                                                                                             ###
###################################################################################################################

$VERSION        =       1.00;
@ISA            =       qw(Exporter);
@EXPORT         =       qw(tbstrains);

sub tbstrains {
	# Switches autoflush for direct printing on.
	$| = 1; 
	# Get parameter and input from front-end.
	my $VAR_dir			=	shift;
	my $POS_OUT			=	shift;
	my $STRAIN_OUT			=	shift;
	my $ref                         =       shift;
	my $refg			=	shift;
	my $micovf			= 	shift;
	my $micovr			=	shift;
	my $mifreq			=	shift;
	my $miphred20			=	shift;
	my $date_string			=	shift;
	my @position_tables		=	@_;
	my $output_mode                 =       3;
	my $genes			=	{};
	my $annotation			=	{};
	my $phylo_positions_homolka	=	{};
	my $phylo_positions_coll	=	{};
	my $phylo_positions_beijing	=	{};
	my $phylo_positions             =       {};
	my $positions_data		=	{};
	my %check_up;
	my $output_file			=	"strain_classification.tab";
	# Check what is already existing in this file.
	if(-f "$STRAIN_OUT/$output_file") {
		open(IN,"$STRAIN_OUT/$output_file") || die "<ERROR>\t",timer(),"\tCan't open $output_file: $!\n";
		<IN>;
		while(<IN>) {
			my $line		=	$_;
			$line           	=~      s/\015?\012?$//;
			my @fields		=	split(/\t/);
			$check_up{$fields[1]."_".$fields[2]."_".$fields[3]."_".$fields[4]}	=	$fields[0];
		}
	}
	# Parse the genomic sequence for determining substitutions.
        print  "<INFO>\t",timer(),"\tParsing $ref...\n";
        my $genome                      =       parse_fasta($VAR_dir,$ref);
        print  "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
	print  "<INFO>\t",timer(),"\tParsing $refg...\n";
        # Parse annotation.
	parse_annotation($VAR_dir,$genes,$annotation,$refg);
        print  "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
	print  "<INFO>\t",timer(),"\t","Start loading classification schemes...\n";
	# Parse lineage classification.
	parse_classification($phylo_positions_homolka,$phylo_positions_coll,$phylo_positions_beijing);
	print  "<INFO>\t",timer(),"\t","Finsihed loading classification schemes!\n";
	print  "<INFO>\t",timer(),"\t","Start combining classification schemes...\n";
	# This hash will contain the data we need.
	foreach my $key(sort { $a cmp $b } keys %$phylo_positions_homolka) {
		$phylo_positions->{$key}->{0}		=	$phylo_positions_homolka->{$key};
	}
	foreach my $key(sort { $a cmp $b } keys %$phylo_positions_coll) {
		$phylo_positions->{$key}->{0}		=	$phylo_positions_coll->{$key};
	}
	foreach my $key(sort { $a cmp $b } keys %$phylo_positions_beijing) {
		$phylo_positions->{$key}->{0}		=	$phylo_positions_beijing->{$key};
	}
	print  "<INFO>\t",timer(),"\t","Finished combining classification schemes!\n";
	print  "<INFO>\t",timer(),"\t","Start writing $output_file...\n";
	unless(-f "$STRAIN_OUT/$output_file") {
		open(OUT,">$STRAIN_OUT/$output_file") || die "<ERROR>\t",timer(),"\tCan't create $output_file: $!\n";
		my $header      =       "Date\tSampleID\tLibraryID\tSource\tRun";
        	$header         .=      "\tHomolka species\tHomolka lineage\tHomolka group\tQuality";
        	$header         .=      "\tColl lineage (branch)\tColl lineage_name (branch)\tColl quality (branch)";
        	$header         .=      "\tColl lineage (easy)\tColl lineage_name (easy)\tColl quality (easy)";
        	$header         .=      "\tBeijing lineage (easy)\tBeijing quality (easy)";
        	$header         .=      "\n";
        	print OUT $header;
		close(OUT);
	}
	open(OUT,">>$STRAIN_OUT/$output_file") || die "<ERROR>\t",timer(),"\tCan't create $output_file: $!\n";
	# Start logic...
	foreach my $file(sort { $a cmp $b } @position_tables) {
    		print   "<INFO>\t",timer(),"\t","Start parsing $file...\n";
		$file =~ /(\S+).gatk_position_table\.tab$/;
		my $id 		=	$1;
		my @sample	=	split(/_/,$id);
		if(exists $check_up{"\'$sample[0]"."_"."\'$sample[1]"."_"."\'$sample[2]"."_"."\'$sample[3]"}) {
			print  "<INFO>\t",timer(),"\t","Skipping $file. Lineage classification already existing!\n";
			next;
		}
		# Parse position table file to hash.
		my $position_table		=	{};
		my $phylo_position_table	= 	{};
		parse_position_table($POS_OUT,$file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
		# Skip every position not included in phylo positions.
		foreach my $pos(sort { by_number() } keys %$phylo_positions) {
			foreach my $index(sort { by_number() } keys %{$phylo_positions->{$pos}}) {
				$phylo_position_table->{$pos}->{$index}		=	$position_table->{$pos}->{$index};
			}
		}
		$position_table 		=	{};
		print   "<INFO>\t",timer(),"\t","Finished parsing $file!\n";
		print   "<INFO>\t",timer(),"\t","Start calling variants for $id!\n";
		my $variants 			=	{};
		my $statistics			=	{};
		call_variants($phylo_position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$output_mode);
		$statistics			=	{};
		$positions_data->{$id} = $variants;
		print   "<INFO>\t",timer(),"\t","Finished calling variants for $id!\n";
		print   "<INFO>\t",timer(),"\t","Start writing lineage result for $id...\n";
		foreach my $id (sort { $a cmp $b } keys %$positions_data) {
			my $lineage		=	"unknown";
			my $species		=	"unknown";
			my $lineage_name	=	"unknown";	
			my $quality_homolka	=	"good";
			my $quality_coll	=	"good";
			my $quality_beijing	=	"good";
			my $IDpositions		=	{};
			foreach my $pos(sort { by_number() } keys %$phylo_positions_homolka) {
				my $allel1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
				my $freq1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];		
				my $count1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
				$quality_homolka	=	"ugly"	unless($allel1 =~ /[AGCTagct]/);
				$quality_homolka	=	"bad"	unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
				$IDpositions->{$pos}	=	$allel1;
			}
			foreach my $pos(sort { by_number() } keys(%$phylo_positions_coll)) {
				my $allel1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
				my $freq1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];
				my $count1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
				$quality_coll		=	"ugly"	unless($allel1 =~ /[AGCTagct]/);
				$quality_coll		=	"bad"	unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
				$IDpositions->{$pos} 	= 	$allel1;
			}	
			foreach my $pos(sort { by_number() } keys(%$phylo_positions_beijing)) {
				my $allel1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
				my $freq1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];
				my $count1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
				$quality_beijing	=	"ugly"	unless($allel1 =~ /[AGCTagct]/);
				$quality_beijing	=	"bad"	unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
				$IDpositions->{$pos} 	= 	$allel1;
			}	
			my $quality			=	0;
			my $outline			= 	"\'$date_string\t\'$sample[0]\t\'$sample[1]\t\'$sample[2]\t\'$sample[3]";
			($species,$lineage_name)	=	specificator_homolka($IDpositions);
			$quality			=	$quality_homolka;
			$lineage			=	translate_homolka2coll($lineage_name);
			$species			=	"\'" . $species;
			$lineage                        =       "\'" . $lineage;
			$lineage_name                   =       "\'" . $lineage_name;
			$quality                        =       "\'" . $quality;
			$outline 			.=	"\t$species\t$lineage\t$lineage_name\t$quality";
			($species,$lineage)		=	specificator_coll_branch($IDpositions);
			$quality			=	$quality_coll;
			$lineage_name			=	translate_coll2homolka($lineage);
			$species                        =       "\'" . $species;
      	          	$lineage                        =       "\'" . $lineage;
       	         	$lineage_name                   =       "\'" . $lineage_name;
                	$quality                        =       "\'" . $quality;
			$outline			.=	"\t$lineage\t$lineage_name\t$quality";
			($species,$lineage)		=	specificator_coll_easy($IDpositions);
			$quality			=	$quality_coll;
			$lineage_name			=	translate_coll2homolka($lineage);
			$species                        =       "\'" . $species;
                	$lineage                        =       "\'" . $lineage;
                	$lineage_name                   =       "\'" . $lineage_name;
                	$quality                        =       "\'" . $quality;
			$outline			.=	"\t$lineage\t$lineage_name\t$quality";
			($species,$lineage)		=	specificator_beijing_easy($IDpositions);
			$quality			=	$quality_beijing;
			$lineage_name			=	$lineage;
                	$lineage                        =       "\'" . $lineage;
                	$quality                        =       "\'" . $quality;
			$outline			.=	"\t$lineage\t$quality";	
			$outline			.=	"\n";
			print OUT $outline;
		}
        	$positions_data              		=       {};
		print   "<INFO>\t",timer(),"\t","Finished writing lineage result for $id!\n";
	}
	close(OUT);
	print  "<INFO>\t",timer(),"\t","Finished writing $output_file!\n";
}


1;
