#!/usr/bin/perl

# tabstop is set to 8.

package TBstrains;

use strict;
use warnings;
use diagnostics;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION        =       1.0.0;
@ISA            =       qw(Exporter);
@EXPORT         =       qw(tbstrains);

sub tbstrains {
	# Get parameter and input from front-end.
	my $logprint			=	shift;
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
	my $all_vars			=	shift;
	$all_vars			=	1; # set to one for fetching all genome positions.
	my $snp_vars			=	shift;
	my $lowfreq_vars		=	shift;
	my @position_tables		=	@_;
	my $genes			=	{};
	my $annotation			=	{};
	my $phylo_positions_homolka	=	{};
	my $phylo_positions_coll	=	{};
	my $phylo_positions_beijing	=	{};
	my $phylo_positions             =       {};
	my $positions_data		=	{};
	my %check_up;
	my $output_file			=	"Strain_Classification.tab";
	# Save already detected strains.
	if(-f "$STRAIN_OUT/$output_file") {
		open(IN,"$STRAIN_OUT/$output_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $output_file: TBstrains line: ", __LINE__ , " \n";
		<IN>;
		while(<IN>) {
			my $line		=	$_;
			$line           	=~      s/\015?\012?$//;
			my @fields		=	split(/\t/);
			$check_up{$fields[1]."_".$fields[2]."_".$fields[3]."_".$fields[4]}	=	$fields[0];
		}
	}
	# Parse the genomic sequence for determining substitutions.
        print $logprint "<INFO>\t",timer(),"\tStart parsing $ref...\n";
        my $genome                      =       parse_fasta($logprint,$VAR_dir,$ref);
        print $logprint "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
	print $logprint "<INFO>\t",timer(),"\tStart parsing $refg...\n";
        # Parse annotation.
	parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
        print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
	print $logprint "<INFO>\t",timer(),"\t","Start loading classification schemes...\n";
	# Parse lineage classification.
	parse_classification($phylo_positions_homolka,$phylo_positions_coll,$phylo_positions_beijing);
	print $logprint "<INFO>\t",timer(),"\t","Finsihed loading classification schemes!\n";
	print $logprint "<INFO>\t",timer(),"\t","Start combining classification schemes...\n";
	# This hash will contain the data we need.
	foreach my $pos (keys %$phylo_positions_homolka) {
		$phylo_positions->{$pos}->{0}		=	$phylo_positions_homolka->{$pos};
	}
	foreach my $pos (keys %$phylo_positions_coll) {
		$phylo_positions->{$pos}->{0}		=	$phylo_positions_coll->{$pos};
	}
	foreach my $pos (keys %$phylo_positions_beijing) {
		$phylo_positions->{$pos}->{0}		=	$phylo_positions_beijing->{$pos};
	}
	print $logprint "<INFO>\t",timer(),"\t","Finished combining classification schemes!\n";
	unless(-f "$STRAIN_OUT/$output_file") {
		print $logprint "<INFO>\t",timer(),"\t","Start writing $output_file...\n";
		open(OUT,">$STRAIN_OUT/$output_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't create $output_file: TBstrains line: ", __LINE__ , " \n";
		my $header      =       "Date\tSampleID\tLibraryID";
        	$header         .=      "\tHomolka species\tHomolka lineage\tHomolka group\tQuality";
        	$header         .=      "\tColl lineage (branch)\tColl lineage_name (branch)\tColl quality (branch)";
        	$header         .=      "\tColl lineage (easy)\tColl lineage_name (easy)\tColl quality (easy)";
        	$header         .=      "\tBeijing lineage (easy)\tBeijing quality (easy)";
        	$header         .=      "\n";
        	print OUT $header;
		close(OUT);
		print $logprint "<INFO>\t",timer(),"\t","Finished writing $output_file!\n";
	}
	open(OUT,">>$STRAIN_OUT/$output_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't create $output_file: TBstrains line: ", __LINE__ , " \n";
	# Start logic...
	foreach my $file (sort { $a cmp $b } @position_tables) {
    		print $logprint "<INFO>\t",timer(),"\t","Start parsing $file...\n";
		$file =~ /(\S+)\.gatk_position_table\.tab$/;
		my $id 		=	$1;
		my @sample	=	split(/_/,$id);
		print "ID: $id\n";
		# Check if strains classification already exists.
		if(exists $check_up{"\'$sample[0]"."_"."\'$sample[1]"}) {
			print $logprint "<INFO>\t",timer(),"\t","Skipping $file. Lineage classification already existing!\n";
			next;
		}
		# Parse position table file to hash.
		my $position_table		=	{};
		my $phylo_position_table	= 	{};
		parse_position_table($logprint,$POS_OUT,$file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
		# Skip every position not included in phylo positions.
		foreach my $pos (keys %$phylo_positions) {
			foreach my $index (keys %{$phylo_positions->{$pos}}) {
				$phylo_position_table->{$pos}->{$index}		=	$position_table->{$pos}->{$index};
			}
		}
		$position_table 		=	{};
		print $logprint "<INFO>\t",timer(),"\t","Finished parsing $file!\n";
		print $logprint "<INFO>\t",timer(),"\t","Start fetching variants for $id...\n";
		my $variants 			=	{};
		my $statistics			=	{};
		call_variants($logprint,$phylo_position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
		$statistics			=	{};
		$phylo_position_table		=	{};
		$positions_data->{$id} = $variants;
		print $logprint "<INFO>\t",timer(),"\t","Finished fetching variants for $id!\n";
		print $logprint "<INFO>\t",timer(),"\t","Start writing lineage result for $id...\n";
		foreach my $id (sort { $a cmp $b } keys %$positions_data) {
			my $lineage		=	"unknown";
			my $species		=	"unknown";
			my $lineage_name	=	"unknown";	
			my $quality_homolka	=	"good";
			my $quality_coll	=	"good";
			my $quality_beijing	=	"good";
			my $IDpositions		=	{};
			foreach my $pos (keys %$phylo_positions_homolka) {
				my $allel1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
				my $freq1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];		
				my $count1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
				$quality_homolka	=	"ugly"	unless($allel1 =~ /[AGCTagct]/);
				$quality_homolka	=	"bad"	unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
				$IDpositions->{$pos}	=	$allel1;
			}
			foreach my $pos (keys %$phylo_positions_coll) {
				my $allel1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
				my $freq1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];
				my $count1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
				$quality_coll		=	"ugly"	unless($allel1 =~ /[AGCTagct]/);
				$quality_coll		=	"bad"	unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
				$IDpositions->{$pos} 	= 	$allel1;
			}	
			foreach my $pos (keys %$phylo_positions_beijing) {
				my $allel1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
				my $freq1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];
				my $count1		=	(split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
				$quality_beijing	=	"ugly"	unless($allel1 =~ /[AGCTagct]/);
				$quality_beijing	=	"bad"	unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
				$IDpositions->{$pos} 	= 	$allel1;
			}	
			my $quality			=	0;
			my $outline			= 	"\'$date_string\t\'$sample[0]\t\'$sample[1]";
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
		$variants				=	{};
		print $logprint "<INFO>\t",timer(),"\t","Finished writing lineage result for $id!\n";
	}
	close(OUT);
	undef(%check_up);
}


1;
