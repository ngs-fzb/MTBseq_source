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

package TBamend;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION	=	1.2.0;
@ISA		=	qw(Exporter);
@EXPORT		=	qw(tbamend);

sub tbamend {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $JOIN_OUT		=	shift;
	my $AMEND_OUT		=	shift;
	my $micovf		=	shift;
	my $micovr		=	shift;
	my $miphred20		=	shift;
	my $mifreq		=	shift;
	my $unambigous		=	shift;
	my $window		=	shift;
	my $categories		=	shift;
	my $resi_list_master	=	shift;
	my $int_regions		=	shift;
	my @join_files		=	@_;
	my $cats		=	{};
        my $variant_infos	=	{};
        my $resi_gene		=	{};
	print $logprint "<INFO>\t",timer(),"\tStart parsing $categories..\n";
	parse_categories($logprint,$categories,$cats);
	print $logprint "<INFO>\t",timer(),"\tFinished parsing $categories!\n";
	print $logprint "<INFO>\t",timer(),"\tStart parsing $resi_list_master and $int_regions..\n";
	parse_variant_infos($logprint,$variant_infos,$resi_gene,$resi_list_master,$int_regions);
	print $logprint "<INFO>\t",timer(),"\tFinished parsing $resi_list_master and $int_regions!\n";
	# Start logic...
	foreach my $join_file (sort { $a cmp $b } @join_files) {
		# At present, process file per line to save RAM.
		my ($amended_file,$phylo_file) = amend_joint_table($logprint,$JOIN_OUT,$AMEND_OUT,$cats,$variant_infos,$resi_gene,$unambigous,$micovf,$micovr,$mifreq,$miphred20,$join_file);
		print $logprint "<INFO>\t",timer(),"\tStart creating $phylo_file with window length $window...\n";
		filter_wlength($logprint,$AMEND_OUT,$phylo_file,$window);
		print $logprint "<INFO>\t",timer(),"\tFinished creating $phylo_file with window length $window!\n";
	}
}


1;
