#!/usr/bin/perl

# tabstop is set to 8.

package TBlist;

use strict;
use warnings;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION	=	1.0.0;
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
