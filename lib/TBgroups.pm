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

package TBgroups;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###                                                                                                             ###
### Description: This package creates groups from a joint SNP analyis. It uses the Joint SNP tables that were   ###
### filtered according to phylogeny SNPs                                                         		###
###                                                                                                             ###
### Input:      _amended_u95_phylo_w12.tab                                                                      ###
### Output:     .matrix .groups                                                                                 ###
###                                                                                                             ###
###################################################################################################################

$VERSION        =       1.10;
@ISA            =       qw(Exporter);
@EXPORT         =       qw(tbgroups);


sub tbgroups {
	# Get parameter and input from front-end.
	my $logprint			=	shift;
	my $AMEND_OUT			=	shift;
	my $GROUPS_OUT			=	shift;
	my $distance			=	shift;
	my @joint_files			=	@_;
	my $strip_ids			=       1;
	# Start logic...
	foreach my $joint_file (sort { $a cmp $b } @joint_files) {
		my $matrix_file				=	"";
		my $group_file				=	"";
		my $pivot_hash				=	{};
		my $position_info			=	{};
		my $distance_matrix			=	{};
		my $strains				=	[];
		$matrix_file				=	strip($joint_file,".matrix");
		$group_file				=	strip($joint_file,("_d" . $distance . ".groups"));
		print $logprint "<INFO>\t",timer(),"\tStart parsing $joint_file...\n";
		parse_amend_table($logprint,$AMEND_OUT,$joint_file,$pivot_hash,$strains,$position_info);
		print $logprint "<INFO>\t",timer(),"\tFinsiehd parsing $joint_file!\n";
 		print $logprint "<INFO>\t",timer(),"\tStart building distance matrix for $joint_file...\n";
		build_matrix($GROUPS_OUT,$matrix_file,$pivot_hash,$distance_matrix,$strains);
		print $logprint "<INFO>\t",timer(),"\tFinished building distance matrix for $joint_file!\n";
		print $logprint "<INFO>\t",timer(),"\tStart calling groups for $joint_file...\n";
		call_groups($GROUPS_OUT,$group_file,$distance,$strip_ids,$strains,$distance_matrix);
		print $logprint "<INFO>\t",timer(),"\tFinished calling groups for $joint_file!\n";
	}
}


1;
