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

package TBlist;

use strict;
use warnings;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION	=	1.2.0;
@ISA		= 	qw(Exporter);
@EXPORT		= 	qw(tblist);

sub tblist {
	# Get parameter and input from front-end.
	my $logprint		=	shift;
	my $VAR_dir		=	shift;
	my $MPILE_OUT		=	shift;
	my $POS_OUT		=	shift;
	my $ref                 =       shift;
        my $mibqual             =       shift;
	my $threads		=	shift;
	my @mpileup_files	=	@_;
	my $ref_hash		=	{};
	if(-f "$VAR_dir/$ref") {
		print $logprint "<INFO>\t",timer(),"\tParsing reference genome $ref...\n";
		# Parse the reference .fasta to fill positions not contained in the .mpileup file (eg. uncalled or uncovered positions).
		parse_reference($logprint,$VAR_dir,$ref,$ref_hash);
		print $logprint "<INFO>\t",timer(),"\tReference genome size (bp): ",scalar(keys %$ref_hash),"!\n";
	}
	else {
		print $logprint "<ERROR>\t",timer(),"\tReference genome, $ref not found: TBlist line 70.\n";
		exit 1;
	}
	# Start logic...
	foreach my $mpileup_file (sort { $a cmp $b } @mpileup_files) {
		my $output	=	$mpileup_file;
		$output		=~	s/.mpileup/_position_table.tab/;
		# Create a position table.
		parse_mpile($logprint,$MPILE_OUT,$POS_OUT,$output,$mpileup_file,$ref_hash,$mibqual,$threads);
	}
	@mpileup_files          =       ();
	$ref_hash		=	{};
}


1;
