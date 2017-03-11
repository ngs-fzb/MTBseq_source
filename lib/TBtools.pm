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

package TBtools;

use strict;
use warnings;
use File::Copy;
use List::Util qw(max);
use Statistics::Basic qw(median mean);
use MCE;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

###################################################################################################################
###                                                                                                             ###
### Description: This package is a repository of functions that are used within the pipeline. subroutines are 	###
### sorted according to their use within this package or within the different pipeline steps.			###
###                                                                                          			###
### Input: Depends on the pipeline step                                                                         ###
### Output: Data structures ready for the pipeline script that is using a subroutine                            ###
###                                                                                                             ###
###################################################################################################################

$VERSION	=	1.10;
@ISA		=	qw(Exporter);
@EXPORT		=	qw(parse_reference parse_fasta parse_mpile parse_position_table parse_annotation parse_variant_infos parse_categories parse_variants parse_amend_table parse_classification print_variants prepare_stats print_joint_table_scaffold print_joint_table print_position_stats help nostep badstep strip by_number timer cat log2 amend_joint_table filter_wlength build_matrix get_seq_len call_variants call_groups translate_homolka2coll translate_coll2homolka specificator_beijing_easy specificator_coll_easy specificator_coll_branch specificator_homolka);



###################################################################
###								###
### 				Local				###
###								###
###################################################################



sub reverse_complement { # Create a reverse complement sequence from input.
	my $dna 		=	shift;
	my $rev_dna		=	reverse($dna);
	$rev_dna		=~ 	tr/ACGTacgt/TGCAtgca/;
	return($rev_dna);
}



sub codon2aa { # Lookup codon table.
	my $logprint		=	shift;
	my $codon		=	shift;
	$codon			=	uc($codon);
	my %g			=	(	'TCA'=>'Ser','TCC'=>'Ser','TCG'=>'Ser', 'TCT'=>'Ser',
						'TTC'=>'Phe','TTT'=>'Phe',
						'TTA'=>'Leu','TTG'=>'Leu',
						'TAC'=>'Tyr','TAT'=>'Tyr',
						'TAA'=>'_',
						'TAG'=>'_',
						'TGC'=>'Cys','TGT'=>'Cys',
						'TGA'=>'_',
						'TGG'=>'Trp',
						'CTA'=>'Leu','CTC'=>'Leu','CTG'=>'Leu','CTT'=>'Leu',
						'CCA'=>'Pro','CCC'=>'Pro','CCG'=>'Pro', 'CCT'=>'Pro',
						'CAC'=>'His','CAT'=>'His',
						'CAA'=>'Gln','CAG'=>'Gln',
						'CGA'=>'Arg','CGC'=>'Arg','CGG'=>'Arg','CGT'=>'Arg',
						'ATA'=>'Ile','ATC'=>'Ile','ATT'=>'Ile','ATG'=>'Met',
						'ACA'=>'Thr','ACC'=>'Thr','ACG'=>'Thr','ACT'=>'Thr',
						'AAC'=>'Asn','AAT'=>'Asn',
						'AAA'=>'Lys','AAG'=>'Lys',
						'AGC'=>'Ser','AGT'=>'Ser',
						'AGA'=>'Arg','AGG'=>'Arg',
						'GTA'=>'Val','GTC'=>'Val','GTG'=>'Val','GTT'=>'Val',
						'GCA'=>'Ala','GCC'=>'Ala','GCG'=>'Ala','GCT'=>'Ala',
						'GAC'=>'Asp','GAT'=>'Asp',
						'GAA'=>'Glu','GAG'=>'Glu',
						'GGA'=>'Gly','GGC'=>'Gly','GGG'=>'Gly','GGT'=>'Gly'
					);
	if(exists $g{$codon}) {
		return $g{$codon};
	}
	else {
		print $logprint "Bad codon \"$codon\" detected!!\n";
		return("Bad");
	}
}



sub ascii_translator { # Translates ascii values into phred scrores.
	my $input 		=	shift;
	my $value 		=	ord($input);
	$value 			=	$value - 33;
	return($value);
}



sub multi_fasta { # Creates a multi fasta file from input.
	my $logprint		=	shift;
	my $AMEND_OUT		=	shift;
	my $phylo_fasta		=	shift;
	my $strain_fasta 	=	shift;
	open(FASTA, ">$AMEND_OUT/$phylo_fasta") || print $logprint "<ERROR>\t",timer(),"\tCan't create $phylo_fasta: $!\n";
	# string for complete FASTA output.
	my $fasta;
	foreach my $header(sort {$a cmp $b } keys %$strain_fasta) {
		my $seq		=	$strain_fasta->{$header};
		$fasta		.=	">$header\n";
		$fasta		.=	"$seq\n";
	}
	print FASTA $fasta;
	close(FASTA);
}



sub fastaids { # Creates an output of plain fasta IDs.
	my $logprint		=	shift;
	my $AMEND_OUT		=	shift;
	my $infasta		=	shift;
	my $outfasta		=	strip($infasta,(".plainIDs.fasta") );
	open(IN,"$AMEND_OUT/$infasta")		||	die print $logprint "<ERROR>\t",timer(),"\tCan't open $infasta: $!\n";
	open(OUT,">$AMEND_OUT/$outfasta")	||	die print $logprint "<ERROR>\t",timer(),"\tCan't create $outfasta: $!\n";
	while(<IN>) {
		my $line		=	$_;
		$line 			=~	s/\015?\012?$//;
		next if($line)		=~	/^#/;
		next unless($line);
		if($line =~ /^>/) {
			my @tags	=	split(/_/, $line);
			$line		=	$tags[0];
        	}
		print OUT $line . "\n";
	}
	close(OUT);
	close(IN);
}



###################################################################
###                                                             ###
###                             Parser                          ###
###                                                             ###
###################################################################



sub parse_reference { # Parse a reference genome and save all posiitions as a unique key.
	my $logprint		=	shift;
	my $VAR_dir		=	shift;
	my $ref			=	shift;
	my $ref_hash		=	shift;
	my $position		=	1;
	open(IN,"$VAR_dir/$ref") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $ref, $!\n";
	while(<IN>) {
		$_		=~ 	s/\015?\012?$//; # This will take care of different line break characters.
		my @line	=	split(//, $_) unless($_ =~ /^>/ );
        	while(my $base = shift(@line) ) {
			$base				= 	lc($base);
			$ref_hash->{$position}->{$base} = 	$base;
			$position++;
        	}
    	}
    close(IN);
}



sub parse_fasta { # Parse a reference genome and save the genome as a whole value.
	my $logprint		=	shift;
	my $VAR_dir		=	shift;
	my $ref			=	shift;
	my $fasta		=	{};
	open(IN,"$VAR_dir/$ref") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $ref: $!\n";
	while(<IN>) {
		chomp;
		my $line	=	$_;
		$line		=~	s/\015?\012?$//; # This will take care of different line break characters
		$fasta->{ref}	.=	$line unless($line =~ /^>/);
    	}
	close(IN);
	return($fasta->{ref});
}



sub parse_mpile { # Parse a .mpileup file.

	##################################################################################################################
	##														##
	## The .mpileup format:												##
	## [1]     [2]     [3]     [4]     [5]     [6]									##
	## ID      pos     ref     cov     bases   quality								##
	##														##
	## .       match forward strand											##
	## ,       match reverse strand											##
	## ACGTN   SNP forward												##
	## acgtn   SNP reverse												##
	## ^       start of a read, next character -33 = mapping quality						##
	## $       end or a read											##
	## *       deleted base (indicated in previous line)								##
	## <       reference skip in read (intron)									##
	## >       reference skip in read (intron)									##
	##														##
	## \+[0-9]+[ACGTNacgtn]+       Insertion between this reference position and the next reference position.	##
	##                             The length of the insertion is given by the integer in the pattern,		##
	##                             followed by the inserted sequence.						##
	## -[0-9]+[ACGTNacgtn]+        Deletion from the reference. The deleted bases will be presented as '*' in	##
	##                             the following lines.								##
	##################################################################################################################

	my $logprint		=	shift;
	my $MPILE_OUT		=	shift;
	my $POS_OUT		=	shift;
	my $output		=	shift;
	my $mpileup_file	=	shift;
	my $ref_hash		=	shift;
	my $mibqual		=	shift;
	my $threads		=	shift;
	my $ref_length		=	scalar(keys(%$ref_hash));
	my %position_table; # This hash is for the final position table.
	my %position_table_insertion; # This hash is for the final position table.
	# Temporary file for parallel processing.
	print $logprint "<INFO>\t",timer(),"\tCreating temporary output file...\n";
	open(OUT,">$POS_OUT/$output\.tmp");
	# This subroutine preserves the output order and takes care of double processing.
	sub preserve_order {
		my %tmp; # Lookup hash.
		my $order_id	=	1; # This order id ensures ordered output.
		# The following sub is a call back function.
		return sub {
			my ($chunk_id,$data)	=	@_;
			$tmp{$chunk_id}		=	$data; # Look up hash for preserving order is set.
			while(1) { # Infinite loop for all of the workers that a coming.
				last unless exists $tmp{$order_id}; # If u havent seen this worker before than call back.
           			print OUT delete $tmp{$order_id++}; # Print worker content and delete directly. Post increment ensures that loop will work with next worker ID
      			}
      			return;
   		};
	}
	# Set multi core environment with $threads available cores.
	# A parallel foreach loop with user requested threads but max 8.
	print $logprint "<INFO>\t",timer(),"\tStart parallel processing...\n";
	my $mce         			=       MCE->new(max_workers => $threads, chunk_size => 1, gather => preserve_order);
	$mce->foreach("$MPILE_OUT/$mpileup_file", sub {
		my ($mce,$chunk_ref,$chunk_id)	=	@_;
		my $line			=	$_;
		my %pile_line; # For saving local line results.
		my %pile_line_insertion; # For saving local line insertions.
		my $index                       =       0; # I need this for non-insertion positions. Because I switched to a three dimensional hash. Now, every insertion will get its pattern.
		$line				=~ 	s/\015?\012?$//; # This will take care of different line break characters.
		my @fields			=	split(/\s+/, $line); # I split according to multiple white spaces because there was an error with tabs.
		# Fields 0 to 2 seem to work every time.
		my $chromosome			=	$fields[0];
		my $position			=	$fields[1];
		my $mpileup_ref			=	$fields[2];
		# But there is inconsistency with 3 to 5.
		my $cov				=	0;
		my $bases			=	0;
		my $qualities			=	0;
		$cov                    	=       $fields[3]	if($fields[3] && $fields[4] && $fields[5]);
		$bases				=	$fields[4]	if($fields[3] && $fields[4] && $fields[5]);
		$qualities			=	$fields[5]	if($fields[3] && $fields[4] && $fields[5]);
		# Starting with removing indication signals from bases string.
		# Remove indicator for read start, and the next character indicating read mapping quality.
		# Remove special characters from $bases.
		$bases				=~	s/\^.{1}//g;
		$bases				=~	s/\$//g;
		$bases				=~	s/<//g;
		$bases				=~	s/>//g;
		# Remove deletion indications.
		while($bases =~ /(-[0-9]+[ACGTNacgtn]+)/g) {
			# Save the length and sequence of the insertion.
			my $deletion		=	$1;
			$deletion		=~	/-([0-9]+)([ACGTNacgtn]+)/g;
			my $deletion_length	=	$1;
			my $deletion_sequence	=	$2;
			my @deletion_sequence	=	split(//, $deletion_sequence);
			# There is a stupid problem if there are non-matching bases (SNPs) directly after a deletion or insertion.
			while(scalar(@deletion_sequence) > $deletion_length) {
				pop(@deletion_sequence);
			}
			$deletion_sequence	=	join("", @deletion_sequence);
			my $deletion_string	=	"-" . $deletion_length . $deletion_sequence;
			# Remove the deletion indicator from $bases string, only remove the first instance.
			$bases			=~	s/$deletion_string//;
			for(my $i = 1;$i <= $deletion_length;$i++) {
				my $base                =   	shift(@deletion_sequence);
				my $deletion_position	=   	$position + $i;
				# We want directional information for gaps, therefore we parse this instead of '*' indicators.
				$pile_line{$deletion_position}{$index}{GAP}++ 	if($base =~ /[ACGTN]/); # The $index is set to 0. Deletions can only occure in the reference.
				$pile_line{$deletion_position}{$index}{gap}++	if($base =~ /[acgtn]/); # The $index is set to 0. Deletions can only occure in the reference.
			}
      	 	}
		# Remove insertion indications.
		while($bases =~ /(\+[0-9]+[ACGTNacgtn]+)/g) {
			# Save the length and sequence of the insertion.
			my $insertion		=   	$1;
			$insertion		=~	/\+([0-9]+)([ACGTNacgtn]+)/g;
			my $insertion_length	=	$1;
			my $insertion_sequence	=	$2;
			my @insertion_sequence	=	split(//, $insertion_sequence);
			# There is a stupid problem if there are non-matching bases directly after an deletion or insertion.
			while(scalar(@insertion_sequence) > $insertion_length) {
	               		pop(@insertion_sequence);
	       	    	}
			$insertion_sequence	= 	join("", @insertion_sequence);
			my $insertion_string	= 	$insertion_length . $insertion_sequence;
			# Remove the insertion indicator from $bases string, only remove the first instance.
			$bases			=~	s/\+$insertion_string//;
			for(my $i = 1;$i <= $insertion_length;$i++) {
	       	        	my $base	= 	shift(@insertion_sequence);
				$pile_line_insertion{$position}{$i}{$base}++		if($base =~ /[ACGTNacgtn]/); # An insertion is outsite the reference. The $index is incremented.
			}
		}
		# Now $bases and $qualities should have equal length.
		my @bases			=	split(//, $bases);
	        my @qualities			=	split(//, $qualities);
		unless(scalar(@bases) == scalar(@qualities)) {
			print $logprint "<WARN>\t",timer(),"\tQualities array and Bases array are not equal in size in genome position $position! Maybe mpileup creation is error proned?\n";
			return(0);
		}
		while(my $base = shift(@bases)) {
			my $quality		=	shift(@qualities);
			my $translated_quality	=	ascii_translator($quality);
			next unless($translated_quality >= $mibqual); # The next is acting on the while loop. Therefore, it is no problem for the entire line and not for MCE.
			$base			=	lc($mpileup_ref)		if($base eq ',');
			$base			=	uc($mpileup_ref)		if($base eq '.');
			$pile_line{$position}{$index}{$base}++				if($base =~ /[ACGTNacgtn]/);
			# Collect quality values.
			if($base =~ /[ACGTNacgtn]/) {
				my $qualbase	=	uc($base);
				my $key_20 	= 	$qualbase . "_qual_20";
			 	$pile_line{$position}{$index}{$key_20}++		if($translated_quality >= 20); # $index is set to 0. We are in the reference.
			}
			if($base =~ /\*/) {
				$pile_line{$position}{$index}{GAP_qual_20}++		if($translated_quality >= 20); # $index is set to 0. We are in the reference.
			}
		}
		my $out				=	""; # Initialize scalar for gathering output in parallel.
		# We did not create a global genome hash because the worker have no permission to manipulate it.
		foreach my $pos(keys %pile_line) {
			my $ref_tmp;
			$ref_tmp		=	$ref_hash->{$pos}->{A}		if(exists $ref_hash->{$pos}->{A}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{C}          if(exists $ref_hash->{$pos}->{C}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{G}          if(exists $ref_hash->{$pos}->{G}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{T}          if(exists $ref_hash->{$pos}->{T}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{N}          if(exists $ref_hash->{$pos}->{N}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{a}          if(exists $ref_hash->{$pos}->{a}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{c}          if(exists $ref_hash->{$pos}->{c}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{g}          if(exists $ref_hash->{$pos}->{g}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{t}          if(exists $ref_hash->{$pos}->{t}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp        	=       $ref_hash->{$pos}->{n}          if(exists $ref_hash->{$pos}->{n}); # We get the information of a base at a certain position from the ref_hash. Subroutine parse_reference was modified for it.
			$ref_tmp		=	uc($ref_tmp);
			my $index		=	0; # The $index is set to 0. We are in the reference.
			# Forward read nucleotides count.
			# Retrive results from local hash.
			my $As			=	0;
			my $Cs			=	0;
			my $Gs			=	0;
			my $Ts			=	0;
			my $Ns			=	0;
			my $GAPs		=	0;
			$As			=	$pile_line{$pos}{$index}{A}	if(exists $pile_line{$pos}{$index}{A});  
			$Cs			=	$pile_line{$pos}{$index}{C}	if(exists $pile_line{$pos}{$index}{C});
			$Gs			=	$pile_line{$pos}{$index}{G}	if(exists $pile_line{$pos}{$index}{G});
			$Ts			=	$pile_line{$pos}{$index}{T}	if(exists $pile_line{$pos}{$index}{T});
			$Ns			=	$pile_line{$pos}{$index}{N}	if(exists $pile_line{$pos}{$index}{N});
			$GAPs			=	$pile_line{$pos}{$index}{GAP}	if(exists $pile_line{$pos}{$index}{GAP});
			# Reverse read nucleotides count.
			my $as			=	0;
			my $cs			=	0;
			my $gs			=	0;
			my $ts			=	0;
			my $ns			=	0;
			my $gaps		=	0;
			$as			=	$pile_line{$pos}{$index}{a}	if(exists $pile_line{$pos}{$index}{a});
			$cs			=	$pile_line{$pos}{$index}{c}	if(exists $pile_line{$pos}{$index}{c});
			$gs			=	$pile_line{$pos}{$index}{g}	if(exists $pile_line{$pos}{$index}{g});
			$ts			=	$pile_line{$pos}{$index}{t}	if(exists $pile_line{$pos}{$index}{t});
			$ns			=	$pile_line{$pos}{$index}{n}	if(exists $pile_line{$pos}{$index}{n});
			$gaps			=	$pile_line{$pos}{$index}{gap}	if(exists $pile_line{$pos}{$index}{gap});
			# Quality Values above Q20.
			my $A_qual_20		=	0;
			my $C_qual_20		=	0;
			my $G_qual_20		=	0;
			my $T_qual_20		=	0;
			my $N_qual_20		=	0;
			my $GAP_qual_20		=	0;
			$A_qual_20		=	$pile_line{$pos}{$index}{A_qual_20}	if(exists $pile_line{$pos}{$index}{A_qual_20});
			$C_qual_20		=	$pile_line{$pos}{$index}{C_qual_20}	if(exists $pile_line{$pos}{$index}{C_qual_20});
			$G_qual_20		=	$pile_line{$pos}{$index}{G_qual_20}	if(exists $pile_line{$pos}{$index}{G_qual_20});
			$T_qual_20		=	$pile_line{$pos}{$index}{T_qual_20}	if(exists $pile_line{$pos}{$index}{T_qual_20});
			$N_qual_20		=	$pile_line{$pos}{$index}{N_qual_20}	if(exists $pile_line{$pos}{$index}{N_qual_20});
			$GAP_qual_20		=	$pile_line{$pos}{$index}{GAP_qual_20}	if(exists $pile_line{$pos}{$index}{GAP_qual_20});
			# Constructing output line.
			my $position_line	=	"$pos\t$index\t$ref_tmp";
			$position_line		.=	"\t$As\t$Cs\t$Gs\t$Ts\t$Ns\t$GAPs";
			$position_line		.=	"\t$as\t$cs\t$gs\t$ts\t$ns\t$gaps";
			$position_line		.=	"\t$A_qual_20\t$C_qual_20\t$G_qual_20\t$T_qual_20\t$N_qual_20\t$GAP_qual_20\n";
			$out			.=	$position_line; # Save the line result to the out scalar that is gathered. If we had deletions in a line, we have more then one position within the hash.
			# If we saw insertions than we can retrive the specific information from this hash.
			if(exists $pile_line_insertion{$pos}) {
				foreach my $insertion_index(keys %{$pile_line_insertion{$pos}}) {
					foreach my $insertion_allel(keys %{$pile_line_insertion{$pos}{$insertion_index}}) {
						$insertion_allel		=       uc($insertion_allel);
						my $Asins			=	0;
						my $Csins			=	0;
						my $Gsins			=	0;
						my $Tsins			=	0;
						my $Nsins			=	0;
						my $GAPsins			=	0;
						my $asins			=	0;
						my $csins			=	0;
						my $gsins			=	0;
						my $tsins			=	0;
						my $nsins			=	0;
						my $gapsins			=	0;
						my $A_qual_20ins		=	0;
						my $C_qual_20ins		=	0;
						my $G_qual_20ins		=	0;
						my $T_qual_20ins		=	0;
						my $N_qual_20ins		=	0;
						my $GAP_qual_20ins		=	0;
						$Asins				=	$pile_line_insertion{$pos}{$insertion_index}{A}	if(exists $pile_line_insertion{$pos}{$insertion_index}{A});
						$Csins				=	$pile_line_insertion{$pos}{$insertion_index}{C}	if(exists $pile_line_insertion{$pos}{$insertion_index}{C});
						$Gsins				=	$pile_line_insertion{$pos}{$insertion_index}{G}	if(exists $pile_line_insertion{$pos}{$insertion_index}{G});
						$Tsins				=	$pile_line_insertion{$pos}{$insertion_index}{T}	if(exists $pile_line_insertion{$pos}{$insertion_index}{T});
						$Nsins				=	$pile_line_insertion{$pos}{$insertion_index}{N}	if(exists $pile_line_insertion{$pos}{$insertion_index}{N});
						$asins				=	$pile_line_insertion{$pos}{$insertion_index}{a}	if(exists $pile_line_insertion{$pos}{$insertion_index}{a});
						$csins				=	$pile_line_insertion{$pos}{$insertion_index}{c}	if(exists $pile_line_insertion{$pos}{$insertion_index}{c});
						$gsins				=	$pile_line_insertion{$pos}{$insertion_index}{g}	if(exists $pile_line_insertion{$pos}{$insertion_index}{g});
						$tsins				=       $pile_line_insertion{$pos}{$insertion_index}{t}	if(exists $pile_line_insertion{$pos}{$insertion_index}{t});
						$nsins				=	$pile_line_insertion{$pos}{$insertion_index}{n}	if(exists $pile_line_insertion{$pos}{$insertion_index}{n});
						# Save also this results to the gather scalar. Now, every insertion is in the right order.
						my $position_line2		=	"$pos\t$insertion_index\t$ref_tmp"; # Unfortunately I cannot save the insertion allel here bcecause joint table module needs here the reference... 
						$position_line2			.=	"\t$Asins\t$Csins\t$Gsins\t$Tsins\t$Nsins\t$GAPsins";
						$position_line2			.=	"\t$asins\t$csins\t$gsins\t$tsins\t$nsins\t$gapsins";
						$position_line2			.=	"\t$A_qual_20ins\t$C_qual_20ins\t$G_qual_20ins\t$T_qual_20ins\t$N_qual_20ins\t$GAP_qual_20ins\n";
						$out				.=	$position_line2;
					}
				}
			}
		}
		@fields		=	(); # Clear @fields to gain RAM.
		undef(%pile_line); # Clear %pileline to gain RAM.
		undef(%pile_line_insertion); # Clear %pile_line_insertion to gain RAM.
		# Go with result directly to preserve order sub.
		MCE->gather($chunk_id,$out);
	});
	$mce->shutdown; # Shuts down parallel processing.
	print $logprint "<INFO>\t",timer(),"\tFinished parallel processing!\n";
	close(OUT);
	# Now we want to fill some gaps that are not in the mpileup file.
	print $logprint "<INFO>\t",timer(),"\tStart loading temporary file into hash structure...\n";
	open(IN,"$POS_OUT/$output\.tmp") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $output\.tmp, $!\n";
        while(<IN>) {
                chomp;
                my @fields	=	split(/\t/, $_);
		my $position	=	shift(@fields);
		my $index	=	shift(@fields);
		my $allel	=	shift(@fields);
		# This is only true for reference and deletions.
		if($index == 0) {
			if(exists $position_table{$position}{$index}{$allel}) {
				my @values=split(/\t/, $position_table{$position}{$index}{$allel});
				for(my $i = 0;$i < scalar(@values);$i++) {
					# This we do because of the parallel processing artifact that we cannot access running childs. Therefore positions with deletions are occuruing more thant one time in the temporary output file.
					$values[$i]	=	$fields[$i] unless($fields[$i] == 0);
				}
				$position_table{$position}{$index}{$allel}		=	join("\t", @values); # more memory efficient than HoA.
			}
			else {
				$position_table{$position}{$index}{$allel}		=	join("\t", @fields); # more memory efficient than HoA.
			}
		}
		else {
			if(exists $position_table_insertion{$position}{$index}{$allel}) {
				my @values=split(/\t/, $position_table_insertion{$position}{$index}{$allel});
				for(my $i = 0;$i < scalar(@values);$i++) {
					$values[$i]     =       $fields[$i] unless($fields[$i] == 0);
				}
				$position_table_insertion{$position}{$index}{$allel}	=       join("\t", @values); # more memory efficient than HoA.
			}
			else {
				$position_table_insertion{$position}{$index}{$allel}    =       join("\t", @fields); # more memory efficient than HoA.
			}
		}
	}
	close(IN);
	print $logprint "<INFO>\t",timer(),"\tFinished loading temporary file into hash structure!\n";
	unlink("$POS_OUT/$output\.tmp") || print $logprint "<WARN>\t",timer(),"\tCan't delete $output\.tmp: no such file!\n";
	print $logprint "<INFO>\t",timer(),"\tStart creating final output file...\n";
	open(OUT,">$POS_OUT/$output");
	my $header		=	"#Pos\tInsindex\tRefBase\tAs\tCs\tGs\tTs\tNs\tGAPs";
        $header			.=	"\tas\tcs\tgs\tts\tns\tgaps";
        $header			.=	"\tAqual_20\tCqual_20\tG_qual_20\tTqual_20\tNqual_20\tGAPqual_20\n";
        print OUT $header;
	# Now we take the ref_hash from the parse_reference sub routine.
	foreach my $pos(sort { by_number() } keys %$ref_hash) {
		my $ref_tmp;
                $ref_tmp        =       $ref_hash->{$pos}->{A}          if(exists $ref_hash->{$pos}->{A});
                $ref_tmp        =       $ref_hash->{$pos}->{C}          if(exists $ref_hash->{$pos}->{C});
                $ref_tmp        =       $ref_hash->{$pos}->{G}          if(exists $ref_hash->{$pos}->{G});
                $ref_tmp        =       $ref_hash->{$pos}->{T}          if(exists $ref_hash->{$pos}->{T});
                $ref_tmp        =       $ref_hash->{$pos}->{N}          if(exists $ref_hash->{$pos}->{N});
                $ref_tmp        =       $ref_hash->{$pos}->{a}          if(exists $ref_hash->{$pos}->{a});
                $ref_tmp        =       $ref_hash->{$pos}->{c}          if(exists $ref_hash->{$pos}->{c});
                $ref_tmp        =       $ref_hash->{$pos}->{g}          if(exists $ref_hash->{$pos}->{g});
                $ref_tmp        =       $ref_hash->{$pos}->{t}          if(exists $ref_hash->{$pos}->{t});
                $ref_tmp        =       $ref_hash->{$pos}->{n}          if(exists $ref_hash->{$pos}->{n});
                $ref_tmp        =       uc($ref_tmp);
		my $index	=	0; # We set index to 0 because we are in the reference.
		my $As		=	0;
		my $Cs		=	0;
		my $Gs		=	0;
		my $Ts		=	0;
		my $Ns		=	0;
		my $GAPs	=	0;
		my $as		=	0;
		my $cs		=	0;
		my $gs		=	0;
		my $ts		=	0;
		my $ns		=	0;
		my $gaps	=	0;
		my $A_qual_20	=	0;
		my $C_qual_20	=	0;
		my $G_qual_20	=	0;
		my $T_qual_20	=	0;
		my $N_qual_20	=	0;
		my $GAP_qual_20	=	0;
		# If we see results in the position_table hash, then we want this information instead.
		if(exists $position_table{$pos}{$index}{$ref_tmp}) {
			print OUT "$pos\t$index\t$ref_tmp\t", delete $position_table{$pos}{$index}{$ref_tmp},"\n";
		}
		# If information was missing, then we print empty values just for completeness.
		else {
			my $position_line	=       "$pos\t$index\t$ref_tmp";
			$position_line          .=      "\t$As\t$Cs\t$Gs\t$Ts\t$Ns\t$GAPs";
			$position_line          .=      "\t$as\t$cs\t$gs\t$ts\t$ns\t$gaps";
			$position_line          .=      "\t$A_qual_20\t$C_qual_20\t$G_qual_20\t$T_qual_20\t$N_qual_20\t$GAP_qual_20\n";
			print OUT $position_line;
		}
		# If we see insertions, then we want to parse another hash instead. Because this hash is outsite of the ref_hash.
		if(exists $position_table_insertion{$pos}) {
			# We have now a true sorted insertion...
			foreach my $insertion_index(sort { by_number() } keys %{$position_table_insertion{$pos}}) {
				# And the insertion allel...
				foreach my $allel(keys %{$position_table_insertion{$pos}{$insertion_index}}) {
					# And we print it after reference position that showed the insertion...
					print OUT "$pos\t$insertion_index\t$allel\t", delete $position_table_insertion{$pos}{$insertion_index}{$allel},"\n";
				}
			}
		}
	}
	undef(%position_table);
	undef(%position_table_insertion);
	close(OUT);
	print $logprint "<INFO>\t",timer(),"\tFinished creating final output file!\n";
	# Runtime approx. 8 minutes on 8 cores.
}



sub parse_position_table { # Parse a position table.
	my $logprint		=	shift;
	my $POS_OUT		=	shift;
	my $position_file	=	shift;
	my $micovf              =       shift;
	my $micovr              =       shift;
	my $miphred20           =       shift;
	my $mifreq              =       shift;
	my $position_table	=	shift;
	my $position_stats	=	shift;
	open(IN,"$POS_OUT/$position_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $position_file: $!\n";
	<IN>;
	while(<IN>) {
		my $line			=	$_;
		$line 				=~	s/\015?\012?$//; # This will take care of different line break characters.
		my @line			=	split(/\t/, $line);
		my $pos				=	shift(@line);
		my $insertion_index		=	shift(@line);
		if($position_stats && $insertion_index == 0) {
			my $ref_tmp		= 	$line[0];
			### Forward read nucleotides count
			my $A			=	$line[1];
			my $C			=	$line[2];
			my $G			=	$line[3];
			my $T			=	$line[4];
			my $N			=	$line[5];
			my $GAP			=	$line[6];
			### Reverse read nucleotides count
			my $a			=	$line[7];
			my $c			=	$line[8];
			my $g			=	$line[9];
			my $t			=	$line[10];
			my $n			=	$line[11];
			my $gap			=	$line[12];
			### Quality Values above Q20
			my $A_qual_20		=	$line[13];
			my $C_qual_20		=	$line[14];
			my $G_qual_20		=	$line[15];
			my $T_qual_20		=	$line[16];
			my $N_qual_20		=	$line[17];
			my $GAP_qual_20		=	$line[18];
			### calculations to fill $statistics_hash
			my $adenosin		=	$A + $a;
			my $cytosin		=	$C + $c;
			my $guanosin		=	$G + $g;
			my $thymin		=	$T + $t;
			my $nucleosin		=	$N + $n;
			my $gaps		=	$GAP + $gap;
			my $coverage		=	$adenosin + $cytosin + $guanosin + $thymin + $nucleosin + $gaps;
			if($coverage == 0) {
				$position_stats->{$pos}->{0}->{nothing}		+=	1;
 	      	 	}
			else {
				my $adenosin_freq	=	$adenosin / $coverage * 100;
				my $cytosin_freq	=	$cytosin / $coverage * 100;
				my $guanosin_freq	=	$guanosin / $coverage * 100;
				my $thymin_freq		=	$thymin / $coverage * 100;
				my $nucleosin_freq	=	$nucleosin / $coverage * 100;
				my $gaps_freq		=	$gaps / $coverage * 100;
				$position_stats->{$pos}->{0}->{any}		+=	1;
				my $unambigous_base_call			= 	0;
				$unambigous_base_call				=	1	if(($A >= $micovf) && ($a >= $micovr) && ($A_qual_20 >= $miphred20) && ($adenosin_freq >= $mifreq));
				$unambigous_base_call				=	1	if(($C >= $micovf) && ($c >= $micovr) && ($C_qual_20 >= $miphred20) && ($cytosin_freq >= $mifreq));
				$unambigous_base_call				=	1	if(($G >= $micovf) && ($g >= $micovr) && ($G_qual_20 >= $miphred20) && ($guanosin_freq >= $mifreq));
				$unambigous_base_call				=	1	if(($T >= $micovf) && ($t >= $micovr) && ($T_qual_20 >= $miphred20) && ($thymin_freq >= $mifreq));
				$unambigous_base_call				=	1	if(($GAP >= $micovf) && ($gap >= $micovr) && ($GAP_qual_20 >= $miphred20) && ($gaps_freq >= $mifreq));
       		    	 	$position_stats->{$pos}->{0}->{unambigous}	+=	1	if($unambigous_base_call);
       		 	}
		}
		$position_table->{$pos}->{$insertion_index}	=	join("\t", @line);
	}
	close(IN);
}



sub parse_annotation { # Parse an annotation file.
	my $logprint		=	shift;
	my $VAR_dir		=	shift;
    	my $genes		=	shift;
	my $annotation		=	shift;
	my $refg		=	shift;
	open(IN,"$VAR_dir/$refg") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $refg: $!\n";
	while(<IN>) {
		chomp;
		my $line			=	$_;
		$line 				=~	s/\015?\012?$//; # This will take care of different line break characters
		next 					if ($line =~ m/^#/);
		my @line			=	split(/\t/, $line);
		my $id                  	=   	$line[0];
		my $name                	=   	$line[1];
		my $start               	=   	$line[2];
		my $stop                	=   	$line[3];
		my $frame               	=   	$line[4];
		my $product         		=   	$line[5];
		my $description         	=   	$line[6];
		my $cogfunction         	=   	$line[7];
		my $cogcats			=   	$line[8];
		my $status_region_name		=   	$line[9];
		my $status_function_name	=   	$line[10];
		my $gene_type			=   	$line[11];
		my $entry			=   	{	'id'			=>	$id,
								'name'			=>	$name,
								'status_region_name'	=>	$status_region_name,
								'status_function_name'	=>	$status_function_name,
								'cogcats'		=>	$cogcats,
								'cogfunction'		=>	$cogfunction,
								'start'			=>	$start,
								'stop'			=>	$stop,
								'product'		=>	$product,
								'description'		=>	$description,
								'frame'			=>	$frame,
								'type'			=>	$gene_type
                                   	 		};
		$genes->{$id}			=	$entry;
		# Reminder: overlapping ORFs on same strand are ignored when using this strategy. This should be changed.
		if($start < $stop) {
			for(my $i = $start;$i <= $stop;$i++) {
				$annotation->{$i}	=	$id;
			}
		}
		else {
			for(my $i = $stop;$i <= $start;$i++) {
				$annotation->{$i}	=	$id;
			}
		}
	}
	close(IN);
}



sub parse_variant_infos { # Parse infos about resistance and phylo SNPs or interesting regions.
	    
	######################################################################
    	##                                                                  ##
    	##  Resi-List-Master.txt columns:                                   ##
    	##                                                                  ##
    	##  0 Variant position genome start                                 ##
    	##  1 Variant position genome stop                                  ##
    	##  2 Var. type                                                     ##
    	##  3 Number                                                        ##
    	##  4 WT base                                                       ##
    	##  5 Var. base                                                     ##
    	##  6 Region                                                        ##
    	##  7 Gene ID                                                       ##
    	##  8 Gene Name                                                     ##
   	##  9 Gene start                                                    ##
    	##  10 Gene stop                                                    ##
    	##  11 Gene length                                                  ##
    	##  12 Dir.                                                         ##
    	##  13 WT AA                                                        ##
    	##  14 Codon nr.                                                    ##
    	##  15 Codon nr. E. coli                                            ##
    	##  16 Var. AA                                                      ##
    	##  17 AA change                                                    ##
    	##  18 Codon change                                                 ##
    	##  19 Variant position gene start                                  ##
    	##  20 Variant position gene stop                                   ##
    	##  21 Antibiotic/Phylogenetic                                      ##
    	##  22 Reference PMID                                               ##
    	##  23 High Confidence SNP                                          ##
    	##                                                                  ##
    	######################################################################
    
	my $logprint		=	shift;
	my $variant_infos	=	shift;
	my $resi_gene		=	shift;
	my $resi_list_master	= 	shift;
	my $int_regions		= 	shift;
	my $genes		= 	shift;
	print $logprint ("<WARN>\t",timer(),"\tNo resistance file $resi_list_master. Will skip resistance annotation.\n")	unless(-f $resi_list_master);
	print $logprint ("<WARN>\t",timer(),"\tNo regions file $int_regions. Will skip resistance annotation.\n")		unless(-f $int_regions);
	open(IN,$resi_list_master) || return;
	# Skip header line.
    	<IN>;
	while(<IN>) {
        	# Delete line break characters.
		$_		=~	s/\015?\012?$//;
        	# Split fields in line, defined by tab seperator.
        	my @line	=	split(/\t/, $_);
        	# Because variants within the $variants hash are upper case, I need to take care about my variant case.
        	$line[5]	=	uc($line[5]);
        	# This step is important because the reference and the variant are in gene orientation within the "Resi-List-Master.txt".
        	# We want to have genome orientation instead.
       		$line[5]	=~	tr/ATGC/TACG/ unless($line[12] eq '+');
		# If we see "SNP".
        	if($line[2] eq 'SNP') {
            		# Then we want to distiguish between resistance- and phylogentic information.
            		if($line[21] =~ /.*phylo.*/) {
				$variant_infos->{$line[0]}->{$line[5]}->{PHYLO}	= 	$line[21];
            		}
            		else {
				# If we are in a resistance position than save the gene name.
                		$resi_gene->{$line[7]}                          = 	$line[8];
                		$variant_infos->{$line[0]}->{$line[5]}->{RESI}	= 	$line[21];
            		}
        	}
        	# If we see "Del".
		if($line[2] eq 'Del') {
            		# Then we want to distiguish between resistance- and phylogentic information.
            		if($line[21] =~ /.*phylo.*/) {
                		for(my $i = ($line[0] + 1);$i <= $line[1];$i++) {
                    			$variant_infos->{$i}->{Del}->{PHYLO}    =	$line[21];
                		}
            		}
            		else {
				# If we are in a resistance determining gene than save the gene name.
                		$resi_gene->{$line[7]}				= 	$line[8];
                		for(my $i = ($line[0] + 1);$i <= $line[1];$i++) {
                    			$variant_infos->{$i}->{Del}->{RESI}	= 	$line[21];
                		}
            		}
        	}
		# If we see an "Ins" entry.
        	if($line[2] eq 'Ins') {
            		# We need the insertion pattern for the entry, therefore we save the information for each specific inserted nuleotide.
            		my @insertion_string		=	split(//, $line[5]);
            		# Then we want to distiguish between resistance- and phylogenetic information.
            		if($line[21] =~ /.*phylo.*/) {
                		for(my $i = 1;$i <= $line[3];$i++) {
                    			$variant_infos->{$line[0]."i".$i}->{$insertion_string[$i]}->{PHYLO}		=	$line[21];
                		}
            		}
            		else {
				# If we are in a resistance determining gene than save the gene name.
                		$resi_gene->{$line[7]}	= 	$line[8];
				my $insertindex		=	1;
                		for(my $i = 1;$i <= $line[3];$i++) {
                    			$variant_infos->{$line[0]."i".$i}->{$insertion_string[$i]}->{RESI}		=	$line[21];
                		}
            		}
        	}
	}
	close(IN);

	######################################################################
	##                                                                  ##
	##  Resi_list_extended_targets_v2.txt columns:                      ##
	##                                                                  ##
	##  0 H37Rv gene                                                    ##
	##  1 start                                                         ##
	##  2 stop                                                          ##
	##  3 type                                                          ##
	##  4 Rv number                                                     ##
	##  5 coding strand                                                 ##
	##  6 associated antibiotic resistance                              ##
	##  7 comment                                                       ##
	##  8 length                                                        ##
	##  9 region analyzed in reference collection                       ##
	##                                                                  ##
	######################################################################

	open(IN,$int_regions) || return;
	# Skip header line.
	<IN>;
    	while(<IN>) {
        	# Delete line break characters.
        	$_          =~  s/\015?\012?$//;
        	# Split fields in line, defined by tab seperator.
        	my @line    =	split(/\t+/, $_);
        	# Save the resistance information entry for all region positions.
		for(my $i = $line[1];$i <= $line[2];$i++) {
            		$variant_infos->{$i}->{TARGET}->{REGION}	=	$line[6];
        	}
    	}
	close(IN);

	# ToDo: Insertions and Deletions are still a problem:
	# Deletion will be annotated if deletion has the following properties:
	#   - deletetion in variant file starts with deletion in Resi-List-Master.txt entry but deletion in variant file is bigger
	#   - deletion of Resi-List-Master.txt entry is nested within a deletion of the variant file
	# Insertions will be annotated if insertions have the following properties:
	#   - they exist exactly like they are in the Resi-List-Master.txt entry but the insertion in the variant file can be also bigger
	#   - they are variable within the inserted sequences and only a few nucleotides are positional identical

}



sub parse_categories { # Parse a category file of essentiell and no-essentiel genes in MTB.
	my $logprint		=	shift;
	my $categories		=	shift;
	my $cats		=	shift;
	print $logprint ("<WARN>\t",timer(),"\tNo categories file $categories. Will skip essential/nonessential annotation.\n") unless(-f $categories);
	open(IN,$categories) || return;
	while(<IN>) {
		my $line	=	$_;
		# Chomp different linebreak characters.
		$line		=~	s/\015?\012?$//;
		next if $line	=~	/^#/;	# This is a useful symbol if you want to uncomment entries within this list.
        	my @line	=	split(/\t/, $line);
        	# Extract ID.
		my $id		=	$line[0];
		my $category	=	$line[1];
		if(exists $cats->{$id}) {
			print $logprint "<WARN>\t",timer(),"\t$id appears twice in input list $categories! Check for inconsistent annotation!\n";
		}
		else {
			$cats->{$id}	=	$category;
		}
	}
	close(IN);
}



sub parse_variants { # Parse a variant file.
	my $logprint		=	shift;
	my $CALL_OUT		=	shift;
	my $file		=	shift;
	my $id			=	shift;
	my $var_positions	=	shift;
	my $strain		=	shift;
	my $micovf		=	shift;
	my $micovr		=	shift;
	my $mifreq		=	shift;
	my $miphred20		=	shift;
	open(IN,"$CALL_OUT/$file") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $file: $!";
	<IN>;
	while(<IN>) {
		chomp;
		my $line			=	$_;
		$line				=~	s/\015?\012?$//; # This will take care of different line break characters.
		my @fields			=	split(/\t/, $line);
		my $pos				=	shift(@fields);
		my $index			=	shift(@fields);
		my $ref_tmp			=	shift(@fields);
		my $type			=	shift(@fields);
		my $allel1			= 	shift(@fields);
		my $covf			= 	shift(@fields);
		my $covr			= 	shift(@fields);
		my $qual20			=	shift(@fields);
		my $freq1			=	shift(@fields);
		my $cov				=	shift(@fields);
		my $subs			=	shift(@fields);
		my $gene			=	shift(@fields);
		my $gene_name			=	shift(@fields);
		my $annotation			=	shift(@fields);
		my $resistance			=	shift(@fields);
		my $phylo			=	shift(@fields);
		my $region			=	shift(@fields);
		my $unambiguous_base_call	=	0;
		# This is for the case of $ouputmode was set to 3 already.
		#if(($covf >= $micovf) && ($covr >= $micovr) && ($freq1 >= $mifreq) && ($qual20 >= $miphred20)) {
		#	if (($allel1 =~ /[ACGTacgt]/) || ($allel1 =~ /GAP/)) {
		#		$unambiguous_base_call = 1;
		#	}
		#}
		#next unless($unambiguous_base_call);
		my @values				=	($type,$allel1,$covf,$covr,$qual20,$freq1,$cov,$subs);
		$var_positions->{$pos}->{$index}	=	$ref_tmp;
		$strain->{$pos.$index.$id}		=	@values;
	}
}



sub parse_amend_table { # Parse an amended joint variant table.
        my $logprint		=	shift;
	my $AMEND_OUT           =       shift;
        my $phylo_file          =       shift;
        my $pivot_hash          =       shift;
        my $strains             =       shift;
        my $position_info       =       shift;
        # Parse phylogeny pivot file.
        open(IN,"$AMEND_OUT/$phylo_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $phylo_file: $!\n";
        # First get the strain header and parse strain names.
        my @strains             =       ();
        my $strain_header       =       <IN>;
        $strain_header          =~      s/\015?\012?$//;
        my @strain_header       =       split(/\t/, $strain_header);
        foreach my $strain(@strain_header) {
                next if $strain         =~      /^#/;
                next unless $strain     =~      /\S+/;
                push(@strains,$strain);
        }
        @$strains               =       @strains;
        # Now get the positions header.
        my $pos_header          =       <IN>;
        while(<IN>) {
                my $line                =       $_;
                # Chomp different linebreak characters.
                $line                   =~      s/\015?\012?$//;
                next unless $line;
                next if $line           =~      /^#/;
                my @line                =       split(/\t/, $line);
                my $pos                 =       shift(@line);
                my $index               =       shift(@line);
                my $ref_tmp             =       shift(@line);
                my $gene                =       shift(@line);
                my $gene_name           =       shift(@line);
                my $product             =       shift(@line);
                my $strains             =       shift(@line);
                my $count_unamb         =       shift(@line);
                my $freq_unamb          =       shift(@line);
                my $Ns                  =       shift(@line);
                my $Us                  =       shift(@line);
                my $PureSNP             =       shift(@line);
                my $main_type           =       shift(@line);
                my $main_allel          =       shift(@line);
                my $category            =       shift(@line);
                my $resistance          =       shift(@line);
                my $phylogeny           =       shift(@line);
                my $region              =       shift(@line);
                # It may be useful to have this information.
                unless(exists($position_info->{$pos}->{$index})) {
                        my @general_values                      =       ($ref_tmp,$gene,$gene_name,$product,$strains,$count_unamb,$Ns,$Us,$PureSNP,$main_type,$main_allel,$category,$resistance,$phylogeny,$region);
                        $position_info->{$pos}->{$index}        =       join("\t", @general_values);
                }
                # Now the line should consist of sets of eight values per strain.
                foreach my $strain(@strains) {
                        my $type                =       shift(@line);
                        my $allel1              =       shift(@line);
                        my $covf                =       shift(@line);
                        my $covr                =       shift(@line);
                        my $qual_20             =       shift(@line);
                        my $freq1               =       shift(@line);
                        my $cov                 =       shift(@line);
                        my $subs                =       shift(@line);
                        my @strain_values       =       ($type,$allel1,$freq1,($covf + $covr),$cov,$subs);
                        $pivot_hash->{$pos}->{$index}->{$strain}        =       join("\t", @strain_values);
                }
        } # At this point, we have finished the current line and variants are recorded.
        close(IN);
}



sub parse_classification {
	my $phylo_positions_homolka			=	shift;
	my $phylo_positions_coll			=	shift;
	my $phylo_positions_beijing			=	shift;
	# homolka2012 classification.
	$phylo_positions_homolka->{648856}		= 	"T";
	$phylo_positions_homolka->{648756}		= 	"C";
	$phylo_positions_homolka->{2053726}		= 	"C";
	$phylo_positions_homolka->{649345}		= 	"C";
	$phylo_positions_homolka->{1128814}		= 	"G";
	$phylo_positions_homolka->{649585}		= 	"G";
	$phylo_positions_homolka->{649601}		= 	"C";
	$phylo_positions_homolka->{157129}		= 	"C";
	$phylo_positions_homolka->{1128825}		= 	"C";
	$phylo_positions_homolka->{1473079}		= 	"G";
	$phylo_positions_homolka->{1674520}		= 	"C";
	$phylo_positions_homolka->{1473094}		= 	"T";
	$phylo_positions_homolka->{1129160}		= 	"C";
	$phylo_positions_homolka->{2053762}		= 	"C";
	$phylo_positions_homolka->{2955957}		= 	"A";
	$phylo_positions_homolka->{648990}		= 	"G";
	$phylo_positions_homolka->{648992}		= 	"C";
	$phylo_positions_homolka->{649067}		= 	"C";
	$phylo_positions_homolka->{157292}		= 	"C";
	$phylo_positions_homolka->{1129165}		= 	"G";
	$phylo_positions_homolka->{2053454}		= 	"G";
	$phylo_positions_homolka->{2956731}		= 	"C";
	$phylo_positions_homolka->{7539}		= 	"A";
	$phylo_positions_homolka->{1816848}		= 	"G";
	$phylo_positions_homolka->{2955233}		=	"C";
	# coll2014 classification.
	$phylo_positions_coll->{615938}			= 	"G";
	$phylo_positions_coll->{4404247}		=	"G";
	$phylo_positions_coll->{3021283}		= 	"G";
	$phylo_positions_coll->{3216553}		=	"G";
	$phylo_positions_coll->{2622402}		=	"G";
	$phylo_positions_coll->{1491275}		=	"G";
	$phylo_positions_coll->{3479545}		= 	"C";
	$phylo_positions_coll->{3470377}		=	"C";
	$phylo_positions_coll->{497491}			=	"G";
	$phylo_positions_coll->{1881090} 		= 	"C";
	$phylo_positions_coll->{2505085}		= 	"G";
	$phylo_positions_coll->{797736}			=	"C";
	$phylo_positions_coll->{4248115}		=	"C";
	$phylo_positions_coll->{3836274}		= 	"G";
	$phylo_positions_coll->{346693} 		=	"G";
	$phylo_positions_coll->{3273107} 		= 	"C";
	$phylo_positions_coll->{1084911} 		= 	"G";
	$phylo_positions_coll->{3722702} 		= 	"G";
	$phylo_positions_coll->{1237818} 		= 	"C";
	$phylo_positions_coll->{2874344} 		= 	"G";
	$phylo_positions_coll->{931123} 		= 	"T";
	$phylo_positions_coll->{62657} 			= 	"G";
	$phylo_positions_coll->{514245} 		= 	"C";
	$phylo_positions_coll->{1850119} 		= 	"C";
	$phylo_positions_coll->{541048} 		= 	"T";
	$phylo_positions_coll->{4229087} 		= 	"C";
	$phylo_positions_coll->{891756} 		= 	"A";
	$phylo_positions_coll->{107794} 		= 	"C";
	$phylo_positions_coll->{2411730} 		= 	"G";
	$phylo_positions_coll->{783601} 		= 	"A";
	$phylo_positions_coll->{1487796} 		= 	"C";
	$phylo_positions_coll->{1455780} 		= 	"T";
	$phylo_positions_coll->{764995} 		= 	"C";
	$phylo_positions_coll->{615614} 		= 	"C";
	$phylo_positions_coll->{4316114} 		= 	"G";
	$phylo_positions_coll->{3388166} 		= 	"C";
	$phylo_positions_coll->{403364} 		= 	"G";
	$phylo_positions_coll->{3977226} 		= 	"G";
	$phylo_positions_coll->{4398141} 		= 	"G";
	$phylo_positions_coll->{1132368} 		= 	"C";
	$phylo_positions_coll->{1502120} 		= 	"C";
	$phylo_positions_coll->{4307886} 		= 	"G";
	$phylo_positions_coll->{4151558} 		= 	"G";
	$phylo_positions_coll->{355181} 		= 	"G";
	$phylo_positions_coll->{2694560} 		= 	"G";
	$phylo_positions_coll->{4246508} 		= 	"G";
	$phylo_positions_coll->{1719757} 		= 	"G";
	$phylo_positions_coll->{3466426} 		= 	"G";
	$phylo_positions_coll->{4260268} 		= 	"G";
	$phylo_positions_coll->{874787} 		= 	"G";
	$phylo_positions_coll->{1501468} 		= 	"G";
	$phylo_positions_coll->{4125058} 		= 	"G";
	$phylo_positions_coll->{3570528} 		= 	"C";
	$phylo_positions_coll->{2875883} 		= 	"C";
	$phylo_positions_coll->{4249732} 		= 	"C";
	$phylo_positions_coll->{3836739} 		= 	"G";
	$phylo_positions_coll->{1759252} 		= 	"G";
	$phylo_positions_coll->{1799921} 		= 	"C";
	$phylo_positions_coll->{1816587} 		= 	"C";
	$phylo_positions_coll->{1137518} 		= 	"G";
	$phylo_positions_coll->{2831482} 		= 	"A";
	$phylo_positions_coll->{1882180}		= 	"C";
	# merker2015 Beijing classification.
	$phylo_positions_beijing->{782634}		=	"A";
	$phylo_positions_beijing->{3061703}		=	"G";
	$phylo_positions_beijing->{2428517}		=	"A";
	$phylo_positions_beijing->{1584762}		=	"G";
	$phylo_positions_beijing->{4238675}		=	"C";
	$phylo_positions_beijing->{46324}		=	"T";
	$phylo_positions_beijing->{15890}		=	"G";
	$phylo_positions_beijing->{4153340}		=	"G";
	$phylo_positions_beijing->{3373397}             =       "G";
}



###################################################################
###                                                             ###
###                             Printer                         ###
###                                                             ###
###################################################################



sub print_variants { # Print a variant file.
	my $CALL_OUT			=	shift;
	my $variants			=	shift;
	my $variant_infos		=	shift;
	my $resi_gene			=	shift;
	my $variants_file		=	shift;
	my $uncovered_file		=	shift;
	my $header			=	"#Pos\tInsindex\tRef\tType\tAllel\tCovFor\tCovRev\tQual20\tFreq\tCov";
	$header				.=	"\tSubst\tGene\tGeneName\tProduct\tResistanceSNP\tPhyloSNP\tInterestingRegion\n";
	my $output_variant		=	$header;
	my $output_uncovered		=	$header;
	open(VAR,">$CALL_OUT/$variants_file");
	open(UNC,">$CALL_OUT/$uncovered_file");
	print VAR $output_variant;
	print UNC $output_uncovered;
	my $gene_name;
	my $resistance;
	my $phylo;
	my $region;
	foreach my $pos(sort { by_number() } keys %$variants) {
		foreach my $insertion_index(sort { by_number() } keys %{$variants->{$pos}}) {
			my @values		=	split(/\t/, $variants->{$pos}->{$insertion_index});
			my $ref_tmp		=	shift(@values);
			my $type		=	shift(@values);
			my $allel1		=	shift(@values);
			my $cov_f		=	shift(@values);
			my $cov_r		=	shift(@values);
			my $qual_20		=	shift(@values);
			my $freq1		=	shift(@values);
			my $cov			=	shift(@values);
			my $subst		=	shift(@values);
			my $gene		=	shift(@values);
			my $gene_name		=	" ";
			$gene_name              =       $resi_gene->{$gene}                             	if(exists $resi_gene->{$gene});
			my $product		=	shift(@values);
			my $resistance		=	" ";
			my $phylo		=	" ";
			my $region		=	" ";
			# SNPs annotation.
			$resistance             =       $variant_infos->{$pos}->{$allel1}->{RESI}		if(exists $variant_infos->{$pos}->{$allel1}->{RESI});
			$phylo                  =       $variant_infos->{$pos}->{$allel1}->{PHYLO}		if(exists $variant_infos->{$pos}->{$allel1}->{PHYLO});
			# Dels annotation.
			$resistance		=	$variant_infos->{$pos}->{Del}->{RESI}			if(exists $variant_infos->{$pos}->{Del}->{RESI});
			$phylo			=	$variant_infos->{$pos}->{Del}->{PHYLO}			if(exists $variant_infos->{$pos}->{Del}->{PHYLO});
			# Ins annotation.
			$resistance		=	$variant_infos->{$pos."i".$insertion_index}->{RESI}	if(exists $variant_infos->{$pos."i".$insertion_index}->{RESI});
			$phylo			=	$variant_infos->{$pos."i".$insertion_index}->{PHYLO}	if(exists $variant_infos->{$pos."i".$insertion_index}->{PHYLO});
			# Region annotation.
			$region                 =       $variant_infos->{$pos}->{TARGET}->{REGION}		if(exists $variant_infos->{$pos}->{TARGET}->{REGION});
			my $line		=	"$pos\t$insertion_index\t$ref_tmp\t$type\t$allel1\t$cov_f\t$cov_r\t$qual_20\t$freq1\t$cov";
			$line			.=	"\t$subst\t$gene\t$gene_name\t$product\t$resistance\t$phylo\t$region\n";
			if($type eq "Unc") {
				print UNC $line;
			}
			else {
				print VAR $line;
			}
		}
	}
	close(VAR);
	close(UNC);
}



sub prepare_stats { # Prints a statistics file.
	my $statistics			=	shift;
	my $ref_genome_size 		=	$statistics->{reference}->{size};
	my $ref_a       		=	$statistics->{reference}->{A};
	my $ref_c       		=	$statistics->{reference}->{C};
	my $ref_g       		=	$statistics->{reference}->{G};
	my $ref_t       		=	$statistics->{reference}->{T};
	my $ref_n       		=	$statistics->{reference}->{N};
	my $ref_gaps    		=	$statistics->{reference}->{gap};
	my $gc_content_ref  		= 	0;
	unless(($ref_a + $ref_t + $ref_g + $ref_c) == 0) {
		$gc_content_ref		=	($ref_g + $ref_c) / ($ref_a + $ref_t + $ref_g + $ref_c) * 100;
	}
	$gc_content_ref			=	sprintf("%.2f",$gc_content_ref);
    	my $any_size    		=	$statistics->{any}->{size};
    	my $any_perc    		=	$any_size / $ref_genome_size;
    	$any_perc       		=	sprintf("%.2f",$any_perc);
	my $any_a       		=	$statistics->{any}->{A};
	my $any_c       		=	$statistics->{any}->{C};
	my $any_g       		=	$statistics->{any}->{G};
    	my $any_t       		=	$statistics->{any}->{T};
    	my $any_n       		=	$statistics->{any}->{N};
    	my $any_gaps    		=	$statistics->{any}->{gap};
    	my $gc_content_any  		= 	0;
    	$gc_content_any 		= 	($any_g + $any_c) / ($any_a + $any_t + $any_g + $any_c) * 100;
    	$gc_content_any 		= 	sprintf("%.2f",$gc_content_any);
    	my @any_coverage 		= 	(0);
	if(exists($statistics->{any}->{coverage})) {
		@any_coverage 		= 	@{$statistics->{any}->{coverage}};
	}
	my $any_mean_coverage   	= 	mean(@any_coverage);
	$any_mean_coverage      	= 	sprintf("%.2f",$any_mean_coverage);
	my $any_median_coverage 	= 	median(@any_coverage);
	my $ein_size    		= 	$statistics->{eindeutig}->{size};
	my $ein_perc    		= 	$ein_size / $ref_genome_size;
	$ein_perc			= 	sprintf("%.2f",$ein_perc);
	my $ein_a       		=	$statistics->{eindeutig}->{A};
	my $ein_c       		=	$statistics->{eindeutig}->{C};
	my $ein_g       		=	$statistics->{eindeutig}->{G};
	my $ein_t       		=	$statistics->{eindeutig}->{T};
	my $ein_n       		=	$statistics->{eindeutig}->{N};
	my $ein_gaps    		=	$statistics->{eindeutig}->{gap};
	my $gc_content_ein      	=	0;
	my $testifnull          	=	($ein_a + $ein_t + $ein_g + $ein_c);
	my $ein_mean_coverage   	=	0;
	my $ein_median_coverage 	=	0;
	if($testifnull > 0) {
		$gc_content_ein		=	($ein_g + $ein_c) / ($ein_a + $ein_t + $ein_g + $ein_c) * 100;
		my @ein_coverage	=	@{$statistics->{eindeutig}->{coverage}};
		$ein_mean_coverage	=	mean(@ein_coverage);
		$ein_median_coverage	=	median(@ein_coverage);
	}
	$gc_content_ein         	=	sprintf("%.2f",$gc_content_ein);
	$ein_mean_coverage		=	sprintf("%.2f",$ein_mean_coverage);
	my $snps			=	$statistics->{SNP};	$snps	=	0	unless($snps);
	my $dels			=	$statistics->{Del};	$dels	=	0	unless($dels);
	my $ins				=	$statistics->{Ins};	$ins	=	0	unless($ins);
	my $unc				=	$statistics->{Unc};	$unc	=	0	unless($unc);
	my $substitutions		=	$statistics->{substitutions};
	$substitutions			=	0						unless($substitutions);
	my $results			=	"'$ref_genome_size\t'$gc_content_ref\t'$any_size ($any_perc)\t'$gc_content_any\t'$any_mean_coverage/$any_median_coverage\t'$ein_size ($ein_perc)\t'$gc_content_ein\t'$ein_mean_coverage/$ein_median_coverage\t";
	$results			=	$results."'$snps\t'$dels\t'$ins\t'$unc\t'$substitutions\n";
	return($results);
}



sub print_joint_table_scaffold { # Prints the scaffold of a joint SNP table.
	my $logprint		=	shift;
	my $JOIN_OUT		=	shift;
	my $join_file		=	shift;
	my $var_positions	=	shift;
	my $annotation		=	shift;
	my $genes		=	shift;
	my @ids			=	@_;
	# Print basic information to output_file.
	open(OUT, ">$JOIN_OUT/$join_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't create $join_file: $!\n";
	# Prepare header for output.
	my $header 		=	join("\t\t\t\t\t\t\t\t", @ids);
	my $second_header 	=	"";
	foreach my $id(@ids) {
		$second_header .= "Type\tAllel\tCovFor\tCovRev\tQual20\tFreq\tCov\tSubst\t";
    	}
	print OUT "#\t\t\t\t\t\t$header\n";
	print OUT "#Position\tInsindex\tRef\tGene\tGeneName\tAnnotation\t$second_header\n";
	# We process input files in turn and print them directly.
	# Print position and annotation info as first step.
	foreach my $pos(sort { by_number() } keys %$var_positions) {
		foreach my $insertion_index(sort { by_number() } keys %{$var_positions->{$pos}}) {
			my $ref_tmp	=	" ";
			my $gene	=	" ";
			my $gene_name	=	" ";
			my $product	=	" ";
			$ref_tmp 	= 	$var_positions->{$pos}->{$insertion_index}		if(exists $var_positions->{$pos}->{$insertion_index});
			$gene       	= 	$annotation->{$pos}					if(exists $annotation->{$pos});
			$gene_name  	= 	$genes->{$gene}->{name}					if(exists $genes->{$gene}->{name});
			$product	= 	$genes->{$gene}->{product}				if(exists $genes->{$gene}->{product});
			print OUT "$pos\t$insertion_index\t$ref_tmp\t$gene\t$gene_name\t$product\n";
		}
    	}
    	close(OUT);
}



sub print_joint_table { # Print a joint SNP table.
	my $logprint		=	shift;
	my $JOIN_OUT		=	shift;
	my $join_file		=	shift;
	my $id			=	shift;
	my $strains		=	shift;
	my $variants		=	shift;
	my $tmp			=	$join_file . ".tmp";
	open(IN, "$JOIN_OUT/$join_file") 	|| die print $logprint "<ERROR>\t",timer(),"\tCan't open $join_file: $!";
	open(OUT, ">$JOIN_OUT/$tmp") 		|| die print $logprint "<ERROR>\t",timer(),"\tCan't create $tmp: $!";
	my $header		=	<IN>;
	my $subheader		=	<IN>;
	print OUT $header;
	print OUT $subheader;
	while(<IN>) {
		my $line	=	$_;
		$line		=~	s/\015?\012?$//; # This will take care of different line break characters.
		my @line	=	split(/\t/, $line);
		my $pos		=	$line[0];
		my $index	=	$line[1];
		my $snpline	=	"";
		my @values	=	split(/\t/, $variants->{$pos}->{$index}) if(exists $variants->{$pos}->{$index});
		my $type	=	$values[1];	$type		=	" "	unless($type);
		my $allel1	=	$values[2];	$allel1		=	" " 	unless($allel1);
		my $covf	=	$values[3];	$covf		=	"0"	unless($covf);
		my $covr	=	$values[4];	$covr		=	"0"	unless($covr);
		my $qual20	=	$values[5];	$qual20		=	"0"	unless($qual20);
		my $freq1	=	$values[6];	$freq1		=	"0"	unless($freq1);
		my $cov		=	$values[7];	$cov		=	"0"	unless($cov);
		my $subs	=	$values[8];	$subs		=	" "	unless($subs);
		$freq1		=	sprintf("%.2f", $freq1) 			unless($freq1 == 0);
		# This marks the information added here not contained in the original variant files
		$allel1		=	lc($allel1)	unless(exists $strains->{$pos.$index.$id});
		$snpline	.=	"\t".$type;
		$snpline	.=	"\t".$allel1;
		$snpline 	.= 	"\t".$covf;
		$snpline	.= 	"\t".$covr;
		$snpline	.= 	"\t".$qual20;
		$snpline	.=	"\t".$freq1;
		$snpline	.=	"\t".$cov;
		$snpline	.=	"\t".$subs;
		print OUT	$line.$snpline."\n";
	}
	close(OUT);
	close(IN);
	move("$JOIN_OUT/$tmp","$JOIN_OUT/$join_file") || die print $logprint "<ERROR>\t",timer(),"\tmove failed: $!\n";
	unlink("$JOIN_OUT/$tmp") || print $logprint "<WARN>\t",timer(),"\tCan't delete $tmp: $!\n";
}



sub print_position_stats {
	my $JOIN_OUT			=	shift;
	my $breadth_file		=	shift;
	my $genome			=	shift;
	my $position_stats		=	shift;
	my @ids				=	@_;
	my $genome_length		=	length($genome);
	my $number_of_datasets		=	scalar(@ids);
	my $number_of_datasets_95	=	int($number_of_datasets / 100 * 95);
	my $number_of_datasets_90	=	int($number_of_datasets / 100 * 90);
	my $number_of_datasets_75	=	int($number_of_datasets / 100 * 75);
	my $positions_covered_100	=	0;
	my $positions_covered_95	=	0;
	my $positions_covered_90	=	0;
	my $positions_covered_75	=	0;
	my $positions_unambigous_100	=	0;
	my $positions_unambigous_95	=	0;
	my $positions_unambigous_90	=	0;
	my $positions_unambigous_75	=	0;
	my $positions_uncovered		=	0;
	foreach my $pos(keys %$position_stats) {
        	my $number_of_ids_with_any_coverage		=	0;
		my $number_of_ids_with_unambiguous_coverage	=	0;
		my $number_of_ids_without_any_coverage		=	0;
		$number_of_ids_with_any_coverage		=	$position_stats->{$pos}->{0}->{any}		if(exists $position_stats->{$pos}->{0}->{any});
		$number_of_ids_with_unambiguous_coverage	=	$position_stats->{$pos}->{0}->{unambigous}	if(exists $position_stats->{$pos}->{0}->{unambigous});
		$number_of_ids_without_any_coverage		=	$position_stats->{$pos}->{0}->{nothing}		if(exists $position_stats->{$pos}->{0}->{nothing});
		$positions_covered_100				+=	1	if($number_of_ids_with_any_coverage == $number_of_datasets);
		$positions_covered_95				+=	1	if($number_of_ids_with_any_coverage >= $number_of_datasets_95);
		$positions_covered_90				+=	1	if($number_of_ids_with_any_coverage >= $number_of_datasets_90);
		$positions_covered_75				+=	1	if($number_of_ids_with_any_coverage >= $number_of_datasets_75);
		$positions_unambigous_100			+=	1	if($number_of_ids_with_unambiguous_coverage == $number_of_datasets);
		$positions_unambigous_95			+=	1	if($number_of_ids_with_unambiguous_coverage >= $number_of_datasets_95);
		$positions_unambigous_90			+=	1	if($number_of_ids_with_unambiguous_coverage >= $number_of_datasets_90);
		$positions_unambigous_75			+=	1	if($number_of_ids_with_unambiguous_coverage >= $number_of_datasets_75);
		$positions_uncovered				+=	1	if($number_of_ids_without_any_coverage == $number_of_datasets);
	}
	open (OUT,">$JOIN_OUT/$breadth_file");
	print OUT "# Number of datasets:\t" . scalar(@ids) . "\n";
	print OUT "# Genome length:\t$genome_length\n";
	print OUT "# Positions covered in 100% of datasets:\t" . $positions_covered_100 . "\n";
	print OUT "# Positions covered in 95% of datasets:\t" . $positions_covered_95 . "\n";
	print OUT "# Positions covered in 90% of datasets:\t" . $positions_covered_90 . "\n";
	print OUT "# Positions covered in 75% of datasets:\t" . $positions_covered_75 . "\n";
	print OUT "# Positions unambigous in 100% of datasets:\t" . $positions_unambigous_100 . "\n";
	print OUT "# Positions unambigous in 95% of datasets:\t" . $positions_unambigous_95 . "\n";
	print OUT "# Positions unambigous in 90% of datasets:\t" . $positions_unambigous_90 . "\n";
	print OUT "# Positions unambigous in 75% of datasets:\t" . $positions_unambigous_75 . "\n";
	print OUT "# Positions uncovered in all datasets:\t" . $positions_uncovered . "\n";
	print OUT "# IDS: \n";
	print OUT join("\n", @ids) . "\n";
	print OUT "\n";
	print OUT "# Position\tInsindex\tUnambigousCoverage\tAnyCoverager\tNoCoverage\n";
	foreach my $pos(sort{ by_number() } keys %$position_stats) {
        my $unambigous_number	=	0;
        my $ambigous_number	=	0;
        my $nothing_number	=	0;
        $unambigous_number      =	$position_stats->{$pos}->{0}->{unambigous}	if(exists $position_stats->{$pos}->{0}->{unambigous});
        $ambigous_number        =	$position_stats->{$pos}->{0}->{any}		if(exists $position_stats->{$pos}->{0}->{any});
        $nothing_number         =	$position_stats->{$pos}->{0}->{nothing}		if(exists $position_stats->{$pos}->{0}->{nothing});
        print OUT "$pos\t0\t$unambigous_number\t$ambigous_number\t$nothing_number\n";
    }
    close(OUT);
}



sub help { # Print a help message.
	print 
	"
	[USAGE]: TBseq [--OPTION PARAMETER]

	
	[DESCRIPTION]: This pipeline generates mappings and calls variants from input samples.

	
	[OPTIONS & PARAMETER]: Please read the README.pdf for detailed information about the parameter!
	
	--step
	<ESSENTIAL> This is an essential option! Choose your pipeline step as a parameter!
		TBfull		Full workflow
		TBreads		Read naming scheme adjustment
		TBbwa		Read mapping
		TBmerge		Merging of duplicate mapping files
		TBrefine	Refinement of mapping(s)
		TBpile		Creation of mpileup file(s)
		TBlist		Creation of position list(s)
		TBvariants	Calling variants
		TBstats         Statisitcs of mapping(s) and variant calling(s)
		TBstrains       Calling lineage from sample(s)
		TBjoin		Joint variant analysis from defined samples
		TBamend		Amending joint variant table(s)
		TBgroups	Detecting groups of samples
	
	--continue
	<OPTIONAL> Ensures that the pipeline continues after selecting a certain pipeline step that is not TBfull.
	
	--samples
	<OPTIONAL> Specifies a column separated text file with sampleID in column 1 and libID in column 2 for pipeline steps after TBstats.

	--project
	<OPTIONAL> Specifiies a project name to identify your joint analysis. Essential for TBamend and TBgroups.

	--resilist
	<OPTIONAL> Specifies a resistance mediating file for resistance prediction. See the README.pdf for file properties.

	--intregions
	<OPTIONAL> Specifies an interesting region files for extended resistance prediction. See the README.pdf for file properties.

	--categories
	<OPTIONAL> Specifies a gene categories file to detect essential and non-essential genes as well as repetitive regions. See the README.pdf for file properties.

	--ref M._tuberculosis_H37Rv_2015-11-13
	<OPTIONAL> Reference genome for mapping. Suported are: \"M._tuberculosis_H37Rv_2015-11-13\", \"M._chimaera_DSM44623_2016-01-28\", \"M._abscessus_CIP-104536T_2014-02-03\" and \"M._fortuitum_CT6_2016-01-08\".
	
	--machine NGS
	<OPTIONAL> Defines the [Source] field within the read naming scheme.

	--run nXXXX
	<OPTIONAL> Defines the [Run] field within the read naming scheme.

	--minbqual 13
        <OPTIONAL> Defines minimum positional mapping quality during variant calling.
	
	--all_vars
	<OPTIONAL> If set, all variant (unambiguous and ambiguous) and invariant sites will be reported.

	--snp_vars
	<OPTIONAL> If set, only unambigous SNPs will be reported. No Insertions nd Deletions will be reported.

	--lowfrew_vars
	<OPTIONAL> If set, alternative low frewuency alleles competing with majority reference alleles will be reported (useful for the detection of subpopulations).

	--mincovf 4
	<OPTIONAL> Defines minimum forward read coverage for a putative variant position.

	--mincovr 4
	<OPTIONAL> Defines minimum reverse read coverage for a putative variant position.

	--minphred20 4
	<OPTIONAL> Defines the minimum number of reads having a phred score above or equal 20 to be considered as a putative variant.

	--minfreq 75
	<OPTIONAL> Defines minimum allele frequency for majority allele.

	--unambig 95
	<OPTIONAL> Defines minimum percentage of samples having unambigous base call in TBamend analysis.

	--window 12
	<OPTIONAL> Defines window for SNP cluster look up. Reduces putative false positives in TBamend.

	--distance 12
	<OPTIONAL> Defines SNP distance for the single linkage clustering in TBgroups.

	--quiet
	<OPTIONAL> Turns off Display logging process.

	--threads 1
	<OPTIONAL> Defines number of CPUs to use. The usage of 8 CPUs is maximum.
	
	
	[EXAMPLES]:

	TBseq --step TBfull
	Default values and execute the whole pipeline.

	TBseq --step TBrefine
	Default values and execute only the \"TBrefine\" step.

	TBseq --step TBbwa --continue
	Default values and execute the \"TBbwa\" module as well as the downstream modules.

	TBseq --step TBfull --threads 8 --machine nextseq --run n0101
	Execute the whole pipeline with 8 threads and supported values.

	TBseq --help
	Print this help message.
	";
	print "\n";
}



sub nostep { # Print an error message when --step is missing.
	print  "\n<ERROR>\t",timer(),"\tNo pipeline step defined. Run: TBseq --help for information!\n\n";
}



sub badstep { # Print an error when --step has a typo.
	my $step	=	shift;
    	print  "\n<ERROR>\t",timer(),"\tBad pipeline step defined! Your parameter [ --step $step ] is no step within the pipeline! Run: start_workflow --help for information!\n\n";
}



#############################################################
###                                                       ###
###                     Helper                            ###
###                                                       ###
#############################################################



sub strip { # Strips the filenme of an input filename.
        my $filename            =       shift;
        my $extension           =       shift;
        $filename               =~      s/\.\w+$//;
        if($extension) {
                $filename       .=      $extension;
        }
        return($filename);
}



sub by_number {
	no warnings;
	$a <=> $b;
}



sub timer { # Sets the time when required.
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	$year	+=	1900;
	$mon	+=	1;
	if(length($mon) == 1 ) {
		$mon		=	"0$mon";
	}
	if(length($mday) == 1 ) {
		$mday		=	"0$mday";
	}
	if(length($hour) == 1) {
		$hour		=	"0$hour";
	}
	if(length($min) ==1 ) {
		$min		=	"0$min";
	}
	if(length($sec) == 1 ) {
		$sec		=	"0$sec";
	}
	my $date_string		=	$year . "-" . $mon . "-" . $mday;
	my $time_string		= 	$hour . ":" . $min . ":" . $sec;
	my $answer_tmp		=	"[" . $date_string . " " . $time_string . "]";
	return($answer_tmp);
}



sub cat { # a simple cat/copy function.
	my $logprint	=	shift;
	my $in		= 	shift;
	my $out		= 	shift;
	open(OUT,">>$out");
	open(IN,"$in") || die print $logprint "<WARN>\t",timer(),"\tCan't open $in, $!\n";
	while(<IN>) {
		print OUT $_;
	}
	close(IN);
	print OUT "\n";
	close(OUT);
}



sub log2 { # returns the log to the base 2.
    my $n 		= 	shift;
    return(log($n)/log(2));
}



sub amend_joint_table { # Amends a joint variant Table.
	my $logprint			=	shift;
	my $JOIN_OUT			=	shift;
	my $AMEND_OUT			=	shift;
	my $cats			=	shift;
	my $variant_infos		=	shift;
	my $resi_gene			=	shift;
	my $unambigous			=	shift;
	my $micovf			=	shift;
	my $micovr			=	shift;
	my $mifreq			=	shift;
	my $miphred20			=	shift;
	my $join_file			=	shift;
	# Create output files.
	my $amended_file		=	strip($join_file,("_amended.tab"));
	my $phylo_file			=	strip($join_file,("_amended_u".$unambigous."_phylo.tab"));
	my $phylo_fasta			=	strip($join_file,("_amended_u".$unambigous."_phylo.fasta"));
	# Open file handle for input file and output files.
	open(IN,"$JOIN_OUT/$join_file")		|| die print $logprint "<ERROR>\t",timer(),"\t Can't open file $join_file: $!\n";
	open(OUT,">$AMEND_OUT/$amended_file")	|| die print $logprint "<ERROR>\t",timer(),"\t Can't create file $amended_file: $!\n";
	open(OUTER,">$AMEND_OUT/$phylo_file")	|| die print $logprint "<ERROR>\t",timer(),"\t Can't create file $phylo_file: $!\n";
	# Get strain header and parse strain names.
	my @strains			=	();
   	my $strain_header		=	<IN>;
	$strain_header			=~	s/\015?\012?$//;
	my @strain_header		=	split(/\t/, $strain_header);
	foreach my $strain(@strain_header) {
		next 			if($strain =~ /^#/);
		next 			unless($strain =~ /\S+/);
		push(@strains,$strain);
	}
   	my $number_of_strains		=	scalar(@strains);
	$strain_header			=	join("\t", @strains);
	# Get the positions header, currently not used.
	my $pos_header			=	<IN>;
   	# Build new header for output file
   	my $first_header		=	join("\t\t\t\t\t\t\t\t",@strains);
   	$first_header			=	"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t".$first_header;
   	my $variant_header		=	"";
   	foreach my $strain(@strains) {
		$variant_header		.=	"\tType\tAllel\tCovFor\tCovRev\tQual20\tFreq\tCov\tSubst"
   	}
	my $second_header		=	"Position\tInsindex\tRef\tGene\tGeneName\tAnnotation";
	$second_header			.=	"\tStrains\tCount_unamb\tPerc_unamb";
	$second_header			.=	"\tNs\tUs\tPureSNP";
	$second_header			.=	"\tMain_allel\tMain_type";
	$second_header			.=	"\tCategory\tResistance\tPhylo\tInterestingRegion";
	$second_header			.=	$variant_header;
   	# After the full output for each position, only the bases and only the amino acids are written in the line.
   	# Add a second tab before strain_header to seperate them from each other.
   	# Maybe better write this to separate files?
	$second_header			.=	"\t\t" . $strain_header . "\t\t" . $strain_header;
	print OUT "#$first_header\n";
	print OUT "#$second_header\n";
	print OUTER "#$first_header\n";
	print OUTER "#$second_header\n";
	# Prepare hash for fasta file.
	my $strain_fasta			=	{};
   	foreach my $strain(@strains) {
       		$strain_fasta->{$strain}	=	"";
   	}
   	# Start processing one line of the pivot file.
   	# We work on each line separately to save RAM.
   	while(<IN>) {
       		my $line			=	$_;
		# Chomp different linebreak characters.
		$line				=~	s/\015?\012?$//;
		next 					unless($line);
		next					if($line =~ /^#/);
		my @line 			=	split(/\t/, $line);
		# Get the four first fields from the joint file.
		my $pos				=	shift(@line);
		my $index			=	shift(@line);
		my $ref_tmp			=	shift(@line);
		my $gene			=	shift(@line);
		my $gene_name			=	shift(@line);
		my $product			=	shift(@line);
		# Collect base information to hash for fasta files.
		my $fasta_hash			=	{};
		# Initialize counters
		my $unambigous_tmp		=	0;
		my $number_Ns			=	0;
		my $number_Us			=	0;
		my $pure_SNP			=	1;
		# Initialize hashes to count base and type.
		my $allel_hash			=	{};
		my $type_hash			=	{};
		# Initialize output lines.
		my $variant_line		=	"";
		my $strain_line			=	"";
		my $codon_line			=	"";
		# Now the line should consist only of sets of eight values for each strain.
		# Collect and sum info for all strains.
		foreach my $strain(@strains) {
			my $type		=	shift(@line);
			my $allel1		=	shift(@line);
			my $covf		=	shift(@line);
			my $covr		=	shift(@line);
			my $qual20		=	shift(@line);
			my $freq1		=	shift(@line);
			my $cov			=	shift(@line);
			my $subs		=	shift(@line);
			# Make sure that frequency is in the right format.
			$freq1			=~	s/,/\./;
			# Count bases and type for majority call.
			unless(($type eq "Unc") || ($allel1 eq " ")) {
				$type_hash->{$type}	+=	1;
				$allel_hash->{$allel1}	+=	1;
			}
			my $unambiguous_base_call	=	0;
			if(($covf >= $micovf) && ($covr >= $micovr) && ($freq1 >= $mifreq) && ($qual20 >= $miphred20)) {
				if(($allel1 =~ /[ACGTacgt]/) || ($allel1 eq "GAP")) {
					$unambiguous_base_call		=	1;
				}
			}
			$unambigous_tmp		+=	1	if($unambiguous_base_call == 1);
			$number_Ns		+=	1	if($allel1 =~ /N/);
			$number_Us		+=	1	if($allel1 =~ /\s/);
			# Concept of pure SNP is important for phylogenetic analysis.
			# - either a SNP or WT (none), no other variant.
			# - base cannot be N or U at any position.
			$pure_SNP		=	0	unless(($type =~ /SNP/) || ($type =~ /none/));
			$pure_SNP		=	0	unless(($number_Ns == 0) && ($number_Us == 0));
			# Collect bases for fasta files
			$fasta_hash->{$strain}->{allel}		=	$allel1		if($pure_SNP == 1);
			$variant_line				.=	"\t$type\t$allel1\t$covf\t$covr\t$qual20\t$freq1\t$cov\t$subs";
			$strain_line				.=	"\t$allel1";
			$codon_line				.=	"\t$subs";
		}
        	# Calculate percentage of strains with unambigous information at this position.
		my $perc_unambigous			=	sprintf("%.2f", $unambigous_tmp / $number_of_strains * 100);
        	# Calculate majority allel and type.
        	my $main_allel				=	" ";
       		my $main_allelcount			=	0;
       		foreach my $allel(keys %$allel_hash) {
           		my $allelcount			=	$allel_hash->{$allel};
           		if($allelcount > $main_allelcount) {
               			$main_allel		=	$allel;
               			$main_allelcount	=	$allelcount;
           		}
        	}
        	my $main_type			=	"Unc";
        	my $main_typecount		=	0;
        	foreach my $type(keys(%$type_hash)) {
          		my $typecount		=	$type_hash->{$type};
			if($typecount > $main_typecount) {
				$main_type	=	$type;
				$main_typecount	=	$typecount;
            		}
        	}
        	# Get further annotation information for this position.
       		my $category			=	" ";
        	my $resistance			=	" ";
        	my $phylo			=	" ";
        	my $region			=	" ";
		my $resistance_gene		=	" ";
		$resistance_gene		=	"yes"							unless(!exists $resi_gene->{$gene});
		$category			=	$cats->{$gene}						unless(!exists $cats->{$gene});
		my $main_allel2			=	uc($main_allel);
        	# SNPs annotation.
		$resistance			=	$variant_infos->{$pos}->{$main_allel2}->{RESI}     	unless(!exists $variant_infos->{$pos}->{$main_allel2}->{RESI});
		$phylo				=	$variant_infos->{$pos}->{$main_allel2}->{PHYLO}    	unless(!exists $variant_infos->{$pos}->{$main_allel2}->{PHYLO});
		# Dels annotation.
		$resistance			=	$variant_infos->{$pos}->{$main_type}->{RESI}     	unless(!exists $variant_infos->{$pos}->{$main_type}->{RESI});
		$phylo				=	$variant_infos->{$pos}->{$main_type}->{PHYLO}     	unless(!exists $variant_infos->{$pos}->{$main_type}->{PHYLO});
		# Ins annotation.
		$resistance			=	$variant_infos->{$pos}->{$main_allel2}->{RESI}     	unless(!exists $variant_infos->{$pos}->{$main_allel2}->{RESI});
		$phylo				=	$variant_infos->{$pos}->{$main_allel2}->{PHYLO}   	unless(!exists $variant_infos->{$pos}->{$main_allel2}->{PHYLO});
		# Interesting region annotation.
		$region				=	$variant_infos->{$pos}->{TARGET}->{REGION}        	unless(!exists $variant_infos->{$pos}->{TARGET}->{REGION});
		# At this point, we can create the output line.
		my $output_line			=	$pos . "\t" . $index ."\t" . $ref_tmp . "\t" . $gene . "\t" . $gene_name . "\t" . $product;
		$output_line			.=	"\t" . $number_of_strains . "\t" . $unambigous_tmp . "\t" . $perc_unambigous;
		$output_line			.=	"\t" . $number_Ns . "\t" . $number_Us . "\t" . $pure_SNP;
		$output_line			.=	"\t" . $main_type . "\t" . $main_allel;
		$output_line			.=	"\t" . $category . "\t" . $resistance . "\t" . $phylo . "\t" . $region;
		$output_line			.=	$variant_line;
		# Add bases and codons separated by another tab.
 		$output_line			.=	"\t" . $strain_line;
		$output_line			.=	"\t" . $codon_line;
		$output_line			.=	"\n";
		print OUT $output_line;
		# For phylogeny output, test for the following criteria.
		# - only true SNPs (variant is either a SNP or WT (none), no other variant, no Ns no Us).
		# - unambigous frequency >= threshold_unambigous.
		# - position not in resistence associated gene.
		# - position not in repetitive genes.
		if(($pure_SNP == 1) && ($perc_unambigous >= $unambigous) && ($category ne "repetitive") && ($resistance_gene ne "yes")) {
			# Print to phylogeny pivot file.
			print OUTER $output_line;
			# Add allel information to fasta hash.
			foreach my $strain(@strains) {
				my $old				=	$strain_fasta->{$strain};
				my $allel			=	$fasta_hash->{$strain}->{allel};
				my $new				=	$old . $allel;
				$strain_fasta->{$strain}	=	$new;
			}
    		}
	} # At this point, we have finished the current line of input and written it to output.
	close(OUT);
	close(OUTER);
	close(IN);
	# Write strain_fasta hash to fasta file and create a plainID fasta file.
	multi_fasta($logprint,$AMEND_OUT,$phylo_fasta,$strain_fasta);
	fastaids($logprint,$AMEND_OUT,$phylo_fasta);
	# Return output file names.
	return($amended_file,$phylo_file);
}



sub filter_wlength { # Fiters for possible fals positive SNPs in a given window.
	my $logprint		=	shift;
	my $AMEND_OUT		=	shift;
	my $phylo_file		=	shift;
	my $window		=	shift;
	my $window_file		=	strip($phylo_file,("_w" . $window . ".tab"));
	my $window_fasta	=	strip($phylo_file,("_w" . $window . ".fasta"));
	my $window_filtered	=	strip($phylo_file,("_w" . $window . "_removed.tab"));
	# Parse phylogeny hash.
	my $pivot_hash		=	{};
	my $strains		=	[];
	my $position_info	=	{};
	parse_amend_table($logprint,$AMEND_OUT,$phylo_file,$pivot_hash,$strains,$position_info);
	my @strains		=	@$strains;
	# Record all positions with SNPs to variant_hash.
	my $variant_hash	=	{};
	foreach my $pos(keys %$pivot_hash) {
		foreach my $index(keys %{$pivot_hash->{$pos}}) {
			foreach my $strain(keys %{$pivot_hash->{$pos}->{$index}}) {
				$variant_hash->{$pos}->{$strain}		=	"present"	if($index == 0);
			}
		}
	}
	# Now use the variant_hash to record variants with a neighbouring SNP within $window.
	# Mark with entry "neighbour" in the window_hash.
	my $window_hash		=	{};
	foreach my $pos(sort { by_number() } keys %$variant_hash) {
		foreach my $strain(sort { $a cmp $b } keys %{$variant_hash->{$pos}}) {
			for(my $i = 1; $i <= $window; $i++) {
				my $index_down	=	$pos - $i;
				my $index_up	=	$pos + $i;
				if((exists($variant_hash->{$index_down}->{$strain})) && (($variant_hash->{$index_down}->{$strain}) eq "present")) {
					$window_hash->{$pos}->{$strain} = "neighbour";
				}
				if((exists($variant_hash->{$index_up}->{$strain})) && (($variant_hash->{$index_up}->{$strain}) eq "present")) {
					$window_hash->{$pos}->{$strain} = "neighbour";
				}
			}
		}
	}
	# Now we know the neighbours and can skip SNPs with neighbours.
	# We parse the phylogeny pivot file and skip everything with a neighbour in $window distance (in the same strain!).
	# At the same time we create the fasta file.
	open(IN,"$AMEND_OUT/$phylo_file")		|| die print $logprint "<ERROR>\t",timer(),"\t Can't open $phylo_file: $!\n";
	open(OUT, ">$AMEND_OUT/$window_file")		|| die print $logprint "<ERROR>\t",timer(),"\t Can't create $window_file: $!\n";
	open(OUTER, ">$AMEND_OUT/$window_filtered")	|| die print $logprint "<ERROR>\t",timer(),"\t Can't create $window_filtered: $!\n";
	my $strain_header_again			=	<IN>;
	my $pos_header_again			=	<IN>;
	print OUT $strain_header_again;
	print OUT $pos_header_again;
	print OUTER $strain_header_again;
	print OUTER $pos_header_again;
	# Prepare to built fasta file of positions sorted for window.
	my $strain_fasta			=	{};
	foreach my $strain(@strains) {
		$strain_fasta->{$strain}	=	"";
	}
	while(<IN>) {
		my $line			=	$_;
		# Save line to tmpline so that we can directly print it.
		my $tmpline			=	$line;
		# Chomp different linebreak characters.
		$line				=~	s/\015?\012?$//;
		next					unless($line);
		next					if($line =~ /^#/);
		my @line			=	split(/\t/, $line);
		my $pos				=	shift(@line);
		my $index			=	shift(@line);
		# Strain specific window filtering.
		my $neighbour			=	0;
		foreach my $strain(@strains) {
			$neighbour		=	1	if((exists($window_hash->{$pos}->{$strain})) && (($window_hash->{$pos}->{$strain}) eq "neighbour"));
		}
		# Now we can print the output if there is no neighbour.
		if($neighbour == 0) {
			print OUT $tmpline;
			# Add allel information to fasta_hash.
			foreach my $strain(@strains) {
				my $old				=	$strain_fasta->{$strain};
				my $allel 			= 	(split(/\t/, $pivot_hash->{$pos}->{$index}->{$strain}))[1];
				my $new 			= 	$old . $allel;
				$strain_fasta->{$strain}	= 	$new;
			}
		}
		if($neighbour > 0) {
			# Print the output to file with removed lines if it was filtered out.
			print OUTER $tmpline;
		}
	} # At this point, we have finished the current line of input and written it to output.
	close(OUTER);
	close(OUT);
	close(IN);
	# Create fasta file and fasta file with clear IDs.
	multi_fasta($logprint,$AMEND_OUT,$window_fasta,$strain_fasta);
	fastaids($logprint,$AMEND_OUT,$window_fasta);
}



sub build_matrix { # Builds a strain matrix.
        my $GROUPS_OUT                  =       shift;
        my $matrix_file                 =       shift;
        my $pivot_hash                  =       shift;
        my $distance_matrix             =       shift;
        my $strains                     =       shift;
	# Initialize matrix.
        foreach my $mainstrain(@$strains) {
                foreach my $strain(@$strains) {
                        $distance_matrix->{$mainstrain}->{$strain} = 0;
                }
        }
        # We calculate only half of the distance matrix, because we do not need more.
        # Remove self-self entries.
        foreach my $mainstrain(keys %$distance_matrix) {
                delete($distance_matrix->{$mainstrain}->{$mainstrain});
        }
        # Remove double entries.
        my @dummy_strains       =       sort(@$strains);
        while(my $mainstrain = shift(@dummy_strains)) {
                foreach my $strain(@dummy_strains) {
                        delete($distance_matrix->{$mainstrain}->{$strain});
                }
        }
        # Fill matrix.
        foreach my $mainstrain(keys %$distance_matrix) {
                foreach my $pos(keys %$pivot_hash) {
                        foreach my $index(keys %{$pivot_hash->{$pos}}) {
                                my $main_allel = uc((split("\t", $pivot_hash->{$pos}->{$index}->{$mainstrain}))[1]);
                                foreach my $strain(keys %{$distance_matrix->{$mainstrain}}) {
                                        my $allel = uc((split("\t", $pivot_hash->{$pos}->{$index}->{$strain}))[1]);
                                        $distance_matrix->{$mainstrain}->{$strain} += 1 unless($main_allel eq $allel);
                                }
                        }
                }
        }
        open(OUT,">$GROUPS_OUT/$matrix_file");
        # Print distance matrix.
        foreach my $mainstrain(sort { $a cmp $b } keys(%$distance_matrix)) {
                print OUT $mainstrain . "\t";
                foreach my $strain(sort { $a cmp $b } keys %{$distance_matrix->{$mainstrain}}) {
                        my $number = $distance_matrix->{$mainstrain}->{$strain};
                        print OUT $number . "\t";
                }
                print OUT "\n";
        }
        print OUT "\n";
        close(OUT);
}



sub get_seq_len {
	my $logprint		=	shift;
	my $W_dir		=	shift;
	my @files 		=	@_;
	my $seq_len;
	open(IN, "gunzip -c $files[0] |") || die print $logprint "<ERROR>\t",timer(),"\t Can't open file $files[0]\n";
		<IN>;
		my $line	=	<IN>;
		chomp($line);
		$seq_len 	=	length($line);
	close(IN);
	$seq_len		.=	"bp";
	return($seq_len);
}



#############################################################
###                                                       ###
###                     Caller                            ###
###                                                       ###
#############################################################



sub call_variants { # Calls variants.
	my $logprint		=	shift;
	my $position_table	=	shift;
	my $variants		=	shift;
	my $statistics		=	shift;
	my $micovf		=	shift;
	my $micovr		=	shift;
	my $miphred20		=	shift;
	my $mifreq		= 	shift;
	my $annotation		=	shift;
	my $genes		=	shift;
	my $ref_genome		=	shift;
	my $all_vars		=	shift;
	my $snp_vars		=	shift;
	my $lowfreq_vars	=	shift;
	foreach my $pos(keys %$position_table) {
		foreach my $insertion_index (keys %{$position_table->{$pos}}) {
			my @values	=	split(/\t/, $position_table->{$pos}->{$insertion_index});
			my $ref_tmp	=	shift(@values);
			my $A		=	shift(@values);
			my $C		=	shift(@values);
			my $G		=	shift(@values);
			my $T		=	shift(@values);
			my $N		=	shift(@values);
			my $GAP		=	shift(@values);
			my $a		=	shift(@values);
			my $c		=	shift(@values);
			my $g		=	shift(@values);
			my $t		=	shift(@values);
			my $n		=	shift(@values);
			my $gap		=	shift(@values);
			my $A_qual_20	=	shift(@values);
			my $C_qual_20	=	shift(@values);
			my $G_qual_20	=	shift(@values);
			my $T_qual_20	=	shift(@values);
			my $N_qual_20	=	shift(@values);
			my $GAP_qual_20	=	shift(@values);
			my $adenosin	=	$A + $a;
			my $cytosin	=	$C + $c;
			my $guanosin	=	$G + $g;
			my $thymin	=	$T + $t;
			my $nucleosin	=	$N + $n;
			my $gaps	=	$GAP + $gap;
			# Coverage (including Ns and gaps).
			my $cov		=	$adenosin + $cytosin + $guanosin + $thymin + $nucleosin + $gaps;
			my $cov_f_all	=	$A + $C + $G + $T + $N + $GAP;
			my $cov_r_all	=	$a + $c + $g + $t + $n + $gap;
			# Call base.
			my $maximum 	=	max($adenosin,$cytosin,$guanosin,$thymin,$nucleosin,$gaps);
			my $allel1	=	" ";
			my $freq1       =	0;
			my $count1	=	0;
			my $qual_20     =	0;
			my $cov_f 	=	0;
			my $cov_r 	=	0;
			my $type        =	"none";
			# In low freq mode, we call the base with highest frequency other than ref base.
			if(($lowfreq_vars == 1) && ($insertion_index == 0) && ($cov > 0)) {
				my @max;
				push(@max,$adenosin)	if(($ref_tmp ne "A") && ($A >= $micovf) && ($a >= $micovr) && ($A_qual_20 >= $miphred20) && ($adenosin / $cov * 100) >= $mifreq);
				push(@max,$cytosin)	if(($ref_tmp ne "C") && ($C >= $micovf) && ($c >= $micovr) && ($C_qual_20 >= $miphred20) && ($cytosin / $cov * 100) >= $mifreq);
				push(@max,$guanosin)	if(($ref_tmp ne "G") && ($G >= $micovf) && ($g >= $micovr) && ($G_qual_20 >= $miphred20) && ($guanosin / $cov * 100) >= $mifreq);
				push(@max,$thymin)	if(($ref_tmp ne "T") && ($T >= $micovf) && ($t >= $micovr) && ($T_qual_20 >= $miphred20) && ($thymin / $cov * 100) >= $mifreq);
				push(@max,$nucleosin)	if(($ref_tmp ne "N") && ($N >= $micovf) && ($n >= $micovr) && ($N_qual_20 >= $miphred20) && ($nucleosin / $cov * 100) >= $mifreq);
				push(@max,$gaps)	if(($ref_tmp ne " ") && ($GAP >= $micovf) && ($gap >= $micovr) && ($GAP_qual_20 >= $miphred20) && ($gaps / $cov * 100) >= $mifreq);
				my $maximum_low		=	max(@max) 	if(scalar(@max) > 0);
				$maximum_low		=	0 		unless($maximum_low);
				# we need to do this in case the reference and a variant have the exact same coverage.
                        	if($maximum_low > 0) {
                                # in case two bases are indicated at the same frequency, it goes A < C < G < T < N < GAP.
                        		if($adenosin == $maximum_low) {
                                        	$allel1		=	"A";
						$count1		=	$adenosin;
                                        	$freq1		=	sprintf("%.2f", $adenosin / $cov * 100);
                                        	$type		=	"SNP" unless($allel1 eq $ref_tmp);
                                        	$qual_20	=	$A_qual_20;
                                        	$cov_f		=	$A;
                                        	$cov_r		=	$a;
                        		}
                        		if($cytosin == $maximum_low) {
                                        	$allel1		=	"C";
						$count1		=	$cytosin;
                                        	$freq1		=	sprintf("%.2f", $cytosin / $cov * 100);
                                        	$type		=	"SNP" unless($allel1 eq $ref_tmp);
						$qual_20	=	$C_qual_20;
                                        	$cov_f		=	$C;
                                       		$cov_r		=	$c;
                        		}
                        		if($guanosin == $maximum_low) {
                                        	$allel1		=	"G";
						$count1		=	$guanosin;
                                        	$freq1		=	sprintf("%.2f", $guanosin / $cov * 100);
                                        	$type		=	"SNP" unless($allel1 eq $ref_tmp);
                                        	$qual_20	=	$G_qual_20;
                                        	$cov_f		=	$G;
                                        	$cov_r		=	$g;
                  			}
                                	if($thymin == $maximum_low) {
                                        	$allel1		=	"T";
						$count1		=	$thymin;
                                        	$freq1		=	sprintf("%.2f", $thymin / $cov * 100);
                                        	$type		=	"SNP" unless($allel1 eq $ref_tmp);
                                        	$qual_20	=	$T_qual_20;
                                        	$cov_f		=	$T;
                                        	$cov_r		=	$t;
                        		}
                        		if($nucleosin == $maximum_low) {
                                        	$allel1		=	"N";
						$count1		=	$nucleosin;
                                        	$freq1		=	sprintf("%.2f", $nucleosin / $cov * 100);
                                        	$type		=	"SNP" unless($allel1 eq $ref_tmp);
                                        	$qual_20	=	$N_qual_20;
                                        	$cov_f		=	$N;
                                        	$cov_r		=	$n;
                        		}
                        		if($gaps == $maximum_low) {
                                        	$allel1		=	"GAP";
						$count1		=	$gaps;
                                        	$freq1		=	sprintf("%.2f", $gaps / $cov * 100);
                                        	$type		=	"Del";
                                        	$qual_20	=	$GAP_qual_20;
                                        	$cov_f		=	$GAP;
                                        	$cov_r		=	$gap;
                        		}
                        	}
			}
			# In case two bases are indicated at the same frequency, it goes A < C < G < T < N < GAP.
			if(($cov > 0) && ($allel1 eq " ")) {
				if($adenosin == $maximum) {
					$allel1		= 	"A";
					$count1         =       $adenosin;
					$freq1		=	sprintf("%.2f", $adenosin / $cov * 100);
					$type		=	"SNP" unless($allel1 eq $ref_tmp);
					$qual_20	=	$A_qual_20;
					$cov_f		=	$A;
					$cov_r		=	$a;
				}
				if($cytosin == $maximum) {
					$allel1		=	"C";
					$count1         =       $cytosin;
					$freq1		=	sprintf("%.2f", $cytosin / $cov * 100);
					$type		=	"SNP" unless($allel1 eq $ref_tmp);
					$qual_20	=	$C_qual_20;
					$cov_f		=	$C;
					$cov_r		=	$c;
				}
				if($guanosin == $maximum) {
					$allel1		=	"G";
					$count1         =       $guanosin;
					$freq1		=	sprintf("%.2f", $guanosin / $cov * 100);
					$type		=	"SNP" unless($allel1 eq $ref_tmp);
					$qual_20	= 	$G_qual_20;
					$cov_f 		= 	$G;
					$cov_r 		= 	$g;
				}
				if($thymin == $maximum) {
					$allel1		=	"T";
					$count1         =       $thymin;
					$freq1		=	sprintf("%.2f", $thymin / $cov * 100);
					$type		=	"SNP" unless($allel1 eq $ref_tmp);
					$qual_20	=	$T_qual_20;
					$cov_f		=	$T;
					$cov_r		=	$t;
				}
				if($nucleosin == $maximum) {
					$allel1 	= 	"N";
					$count1         =       $nucleosin;
					$freq1 		= 	sprintf("%.2f", $nucleosin / $cov * 100);
					$type 		= 	"SNP" unless ($allel1 eq $ref_tmp);
					$qual_20 	= 	$N_qual_20;
					$cov_f		=	$N;
					$cov_r		=	$n;
				}
				if($gaps == $maximum) {
					$allel1 	= 	"GAP";
					$count1         =       $gaps;
					$freq1 		=	sprintf("%.2f", $gaps / $cov * 100);
					$type 		= 	"Del";
					$qual_20 	= 	$GAP_qual_20;
					$cov_f		=	$GAP;
					$cov_r		=	$gap;
				}
			}
			# Uncovered positions.
			unless($cov > 0) {
				$type		=	"Unc";
				$allel1		=	" ";
			}
			# Insertions need special treatment.
			if($insertion_index != 0) {
				$type 			=	"Ins";
				my @values_prev		=	split(/\t/, $position_table->{$pos}->{"0"});
				my $ref_tmp_prev      	=       shift(@values_prev);
				my $A_prev           	=       shift(@values_prev);
                        	my $C_prev           	=       shift(@values_prev);
                        	my $G_prev           	=       shift(@values_prev);
                        	my $T_prev           	=       shift(@values_prev);
                        	my $N_prev           	=       shift(@values_prev);
                        	my $GAP_prev         	=       shift(@values_prev);
                        	my $a_prev           	=       shift(@values_prev);
                        	my $c_prev           	=       shift(@values_prev);
                        	my $g_prev           	=       shift(@values_prev);
                        	my $t_prev           	=       shift(@values_prev);
                        	my $n_prev           	=       shift(@values_prev);
                        	my $gap_prev         	=       shift(@values_prev);
				my $A_qual_20_prev	=       shift(@values_prev);
                        	my $C_qual_20_prev   	=       shift(@values_prev);
                        	my $G_qual_20_prev   	=       shift(@values_prev);
                        	my $T_qual_20_prev   	=       shift(@values_prev);
                        	my $N_qual_20_prev   	=       shift(@values_prev);
                        	my $GAP_qual_20_prev 	=       shift(@values_prev);
                        	my $adenosin_prev    	=       $A_prev + $a_prev;
                        	my $cytosin_prev     	=       $C_prev + $c_prev;
                        	my $guanosin_prev    	=       $G_prev + $g_prev;
                        	my $thymin_prev      	=       $T_prev + $t_prev;
                        	my $nucleosin_prev   	=       $N_prev + $n_prev;
                        	my $gaps_prev        	=       $GAP_prev + $gap_prev;
				my $cov_prev		=       $adenosin_prev + $cytosin_prev + $guanosin_prev + $thymin_prev + $nucleosin_prev + $gaps_prev;
				$qual_20		=	$miphred20;
				if($allel1 eq 'A') {
					$freq1          =       sprintf("%.2f", $adenosin / $cov_prev * 100);
				}
				if($allel1 eq 'C') {
					$freq1          =       sprintf("%.2f", $cytosin / $cov_prev * 100);
				}
				if($allel1 eq 'G') {
					$freq1          =       sprintf("%.2f", $guanosin / $cov_prev * 100);
				}
				if($allel1 eq 'T') {
					$freq1          =       sprintf("%.2f", $thymin / $cov_prev * 100);
				}
				if($allel1 eq 'N') {
					$freq1          =       sprintf("%.2f", $nucleosin / $cov_prev * 100);
				}
			}
			my $unambiguous_base_call		=	0;
			# Strict filtering - every threshold must agree.
			if(($cov_f >= $micovf) && ($cov_r >= $micovr) && ($freq1 >= $mifreq) && ($qual_20 >= $miphred20)) {
				if(($allel1 =~ /[ACGTacgt]/) || ($allel1 eq "GAP")) {
					$unambiguous_base_call	=	1;
				}
			}
			# Statistics part, fill $statistitcs.
			unless($type eq "Ins") {
				# Statistics for reference genome.
				$statistics->{reference}->{size}	+=	1;
				$statistics->{reference}->{$ref_tmp}	+=	1;
				# Statistics covered by any read.
				if($cov > 0) {
					$statistics->{any}->{size}	+=	1;
					$statistics->{any}->{$allel1}	+=	1;
					push(@{$statistics->{any}->{coverage}}, $cov);
					# Statistics unambigous coverage.
					if($unambiguous_base_call == 1) {
						# N is not unambigous.
						if($allel1 ne "N") {
							$statistics->{eindeutig}->{size}	+=	1;
							$statistics->{eindeutig}->{$allel1}	+=	1;
							push(@{$statistics->{eindeutig}->{coverage}}, $cov);
						}
					}
            			}
        		}
			# Counting called variants.
			$statistics->{$type}	+=	1	if($type eq "Unc");
			$statistics->{$type}	+=	1	if(($unambiguous_base_call) && ($type ne "none"));
			if($all_vars != 1) {
				# Skip all positions not different from reference.
				next 		if(($allel1 eq $ref_tmp) && !($type eq "Ins"));
				# Skip all ambigous positions
				next		unless(($unambiguous_base_call == 1) || ($type eq "Unc"));
				next 		if(($snp_vars == 1) && !($type eq "SNP"));
			}
			# Gene specific part.
			my $gene            	=	" ";
			my $gene_type       	=	" ";
			my $product         	=	" ";
			my $substitution    	=	" ";
			my $real_substitution 	=	0;
			$gene			=	$annotation->{$pos}		if(exists $annotation->{$pos});
			$gene_type		=	$genes->{$gene}->{type}		if(exists $genes->{$gene}->{type});
			$product		=	$genes->{$gene}->{product}	if(exists $genes->{$gene}->{product});
			# Calculate substitution if variant is a SNP and in a gene.
			if(($type eq "SNP") && ($gene_type eq "CDS")) {
				my $start		=	$genes->{$gene}->{start};
				my $stop		=	$genes->{$gene}->{stop};
				# Genes in forward direction.
				if($start < $stop) {
					# Length of the gene and pos of SNP in the gene.
					my $length	=	$stop - $start + 1;
					my $gene_pos	=	$pos - $start + 1;
					my $tmp		=	($gene_pos + 2) / 3;
					my $rank	=	int($tmp);
					# Fetch the dna sequence from genome (genome starts with 0).
					my $dna		=	substr($ref_genome, ($start - 1), $length);
					my @codons;
					for(my $i = 0;$i < (length($dna) - 2);$i += 3) {
						my $codon=substr($dna, $i, 3);
						push(@codons,$codon);
					}
					# The codon with the snp is the n-th codon, but array starts at 0.
					my $codon 			= 	$codons[($rank - 1)];
					my $aa 				=	codon2aa($logprint,$codon);
					# SNP dna is like dna with one base exchange.
					my $snp_dna			= 	$dna;
					substr($snp_dna,$gene_pos-1,1)	=	$allel1;
					my @snp_codons;
					for(my $i = 0;$i < (length($snp_dna) - 2); $i += 3) {
						my $codon		=	substr($snp_dna, $i, 3);
						push(@snp_codons,$codon);
					}
					my $snp_codon			=	$snp_codons[($rank - 1)];
					my $snp_aa			=	"bad";
					$snp_aa				=	codon2aa($logprint,$snp_codon) unless($snp_codon =~ /Nn/);
					$real_substitution		=	1 unless($aa eq $snp_aa);
					# Print substitution string even without actual substitution.
					my $aa_tmp_string		=	$aa . $rank . $snp_aa;
					my $codon_tmp_string		= 	"(" . $codon . "/" . $snp_codon . ")";
					$substitution			=	($aa_tmp_string . " " . $codon_tmp_string);
				}
				# Genes in reverse direction.
				if($start > $stop) {
					# Length of the gene and pos of SNP in the gene.
					my $length			=	$start - $stop + 1;
					my $gene_pos			=	abs($pos - $start) + 1;
					my $tmp				=	($gene_pos + 2) / 3;
					my $rank			=	int($tmp);
					my $dna				=	substr($ref_genome, ($stop-1), $length);
					my $rev_dna			=	reverse_complement($dna);
					my $rev_allel1			=	reverse_complement($allel1);
					my $rev_ref			=	reverse_complement($ref_tmp);
					my @codons;
					for(my $i = 0;$i < (length($rev_dna) - 2);$i += 3) {
						my $codon		=	substr($rev_dna, $i, 3);
						push(@codons, $codon);
 	      	         		}
					my $codon				=	$codons[($rank - 1)];
					my $aa					=	codon2aa($logprint,$codon);
					my $snp_dna				=	$rev_dna;
					substr($snp_dna, $gene_pos - 1, 1)	=	$rev_allel1;
					my @snp_codons;
					for(my $i = 0;$i < (length($snp_dna) - 2);$i += 3) {
						my $codon=substr($snp_dna, $i, 3);
						push(@snp_codons,$codon);
                			}
					my $snp_codon			=	$snp_codons[($rank-1)];
					my $snp_aa			=	"bad";
                			$snp_aa				=	codon2aa($logprint,$snp_codon) unless($snp_codon =~ /Nn/);
                			$real_substitution		=	1 unless($aa eq $snp_aa);
                			# Print substitution string even without actual substitution.
                			my $aa_tmp_string		=	$aa . $rank . $snp_aa;
                			my $codon_tmp_string		=	"(" . $codon . "/" . $snp_codon . ")";
                			$substitution			=	($aa_tmp_string . " " . $codon_tmp_string);
				}
			}
			# Statistics part, count substitutions.
			$statistics->{substitutions}		+=	1 if($real_substitution > 0);
			# Fill information into $variants hash.
			@values					=	($ref_tmp,$type,$allel1,$cov_f,$cov_r,$qual_20,$freq1,$cov,$substitution,$gene,$product,$count1);
			$variants->{$pos}->{$insertion_index}	=	join("\t", @values);
		}
	}
}



sub call_groups { # Calls groups from a joint strain analysis.
        my $GROUPS_OUT          =       shift;
        my $group_file          =       shift;
        my $distance            =       shift;
        my $strip_ids           =       shift;
        my $strains             =       shift;
        my $distance_matrix     =       shift;
        my $group_number        =       1;
        my $group_changed       =       1;
        my $max_distance        =       0;
        my $strain_affiliation  =       {};
        my $temp_groups         =       {};
        # The easy approach to call groups, just collect incrementally.
        # Stop if no group changes.
        # Initialize self entry.
        foreach my $mainstrain(keys %$distance_matrix) {
                $strain_affiliation->{$mainstrain} = "ungrouped";
        }
        my $min_distance        =       $distance + 1;
        # We need this to make sure that all strains get the right group .
        # If two groups grow together.
        for(my $i = 0;$group_changed == 1;$i++) {
                $min_distance           =       $distance + 1;
                $group_changed          =       0;
                my $mainstrain_save;
                my $strain_save;
                # Identify pair with smallest distance, take first smallest distance.
                foreach my $mainstrain(keys %$distance_matrix) {
                        foreach my $strain(keys %{$distance_matrix->{$mainstrain}}) {
                                my $number      =       $distance_matrix->{$mainstrain}->{$strain};
                                if($number < $min_distance) {
                                        $min_distance           =       $number;
                                        $mainstrain_save        =       $mainstrain;
                                        $strain_save            =       $strain;
                                }
                        }
                }
                my $group = "temp";
                if($min_distance <= $distance) {
                        $group_changed          =       1;
                        my $mainstrain_group    =       $strain_affiliation->{$mainstrain_save};
                        my $strain_group        =       $strain_affiliation->{$strain_save};
                        if($mainstrain_group =~ /group_\d+/) {
                                $group          =       $mainstrain_group;
                        }
                        if($strain_group =~ /group_\d+/) {
                                $group          =       $strain_group;
                        }
                        if(($mainstrain_group =~ /group_\d+/) && ($strain_group =~ /group_\d+/)) {
                                my $mainstrain_group_number     =       $mainstrain_group       =~      /group_(\d+)/;
                                my $strain_group_number         =       $strain_group           =~      /group_(\d+)/;
                                $group                          =       "group_" . $mainstrain_group_number;
                                $group                          =       "group_" . $strain_group_number if($strain_group_number < $mainstrain_group_number);
                        }
                        unless($group =~ /group_/) {
                                $group                          =       "group_" . "$group_number";
                                $group_number++;
                        }
                        $strain_affiliation->{$mainstrain_save} =       $group;
                        $strain_affiliation->{$strain_save}     =       $group;
                        # We need this to make sure that all strains get the right group.
                        # If two groups grow together.
                        if(($mainstrain_group =~ /group_\d+/) && !($mainstrain_group eq $group)) {
                                my @temp_strains                =       keys(%{$temp_groups->{$mainstrain_group}});
                                foreach my $temp_strain(@temp_strains) {
                                        $strain_affiliation->{$temp_strain}     =       $group;
                                        $temp_groups->{$group}->{$temp_strain}  =       1;
                                }
                                delete($temp_groups->{$mainstrain_group});
                        }
                        if(($strain_group =~ /group_\d+/) && !($strain_group eq $group)) {
                                my @temp_strains        =       keys(%{$temp_groups->{$strain_group}});
                                foreach my $temp_strain(@temp_strains) {
                                        $strain_affiliation->{$temp_strain}     =       $group;
                                        $temp_groups->{$group}->{$temp_strain}  =       1;
                                }
                                delete($temp_groups->{$strain_group});
                        }
                        $temp_groups->{$group}->{$mainstrain_save}      =       1;
                        $temp_groups->{$group}->{$strain_save}          =       1;
                        # We need to do this otherwise we will never get out of the for loop.
                        # Since we work with only half the matrix, we would actually only need one of these.
                        delete($distance_matrix->{$mainstrain_save}->{$strain_save});
                        delete($distance_matrix->{$strain_save}->{$mainstrain_save});
                }
        }
        my $groups              =       {};
        my $groups_renamed      =       {};
        foreach my $strain(@$strains) {
                my $group       =       $strain_affiliation->{$strain};
                if($strip_ids == 1) {
                        my @strain_ids  =       split("_",$strain);
                        $strain         =       $strain_ids[0];
                }
                push(@{$groups->{$group}},$strain);
        }
        # Rename groups to have continous numbering.
        my @groups              =       sort(keys %$groups);
        my $number_of_groups    =       scalar(@groups);
        for (my $i = 1;$i <= $number_of_groups;$i++) {
                my $new_group_name                      =       "group_" . "$i";
                my $old_group_name                      =       shift(@groups);
                $new_group_name                         =       $old_group_name if($old_group_name eq "ungrouped");
                $groups_renamed->{$new_group_name}      =       $groups->{$old_group_name};
        }
        $groups         =       $groups_renamed;
        open(OUT,">$GROUPS_OUT/$group_file");
        print OUT "### Output as groups:\n";
        foreach my $group(sort { $a cmp $b } keys %$groups) {
                my @strains     =       @{$groups->{$group}};
                print OUT "> $group\t" . scalar(@strains) . "\n" . join(",", sort(@strains)) . "\n";
        }
        print OUT "\n\n";
        print OUT "### Output as lists:\n";
        foreach my $group(sort { $a cmp $b } keys %$groups) {
                my @strains = @{$groups->{$group}};
                foreach my $strain (sort(@strains)) {
                        print OUT $strain . "\t" . $group . "\n";
                }
        }
        close(OUT);
}



#############################################################
###                                                       ###
###                     Translator                        ###
###                                                       ###
#############################################################



sub translate_homolka2coll {
	my $homolka_lineage			=	shift;
	my $coll_lineage                        =	'unkown';
	my $homolka_hash			=	{};
	$homolka_hash->{'EAI'}			=	'1';
	$homolka_hash->{'EAI Manila'}		=	'1.2.1';
	$homolka_hash->{'Beijing'}		=	'2';
	$homolka_hash->{'Delhi-CAS'}		=	'3';
	$homolka_hash->{'Haarlem'}		=	'4.1.2.1';
	$homolka_hash->{'Ural'}			=	'4.2.1';
	$homolka_hash->{'TUR'}			=	'4.2.2.1';
	$homolka_hash->{'LAM'}			=	'4.3';
	$homolka_hash->{'S-type'}		=	'4.4.1.1';
	$homolka_hash->{'Uganda'}		=	'4.6.1';
	$homolka_hash->{'Cameroon'}		=	'4.6.2.2';
	$homolka_hash->{'Ghana'}		=	'unknown';
	$homolka_hash->{'West African 1b'}	=	'5';
	$homolka_hash->{'West African 1a'}	=	'5';
	$homolka_hash->{'West African 2'}	=	'6';
	$homolka_hash->{'M. caprae'}		=	'unknown';
	$homolka_hash->{'M. canettii'}		=	'unknown';
	$homolka_hash->{'M. microti'}		=	'unknown';
	$homolka_hash->{'M. pinnipedii'}	=	'unknown';
	$homolka_hash->{'M. bovis'}		=	'BOV';
	$homolka_hash->{'MTBC Clade 2'}		=	'unknown';
	$homolka_hash->{'Clade 1'}		=	'unknown';
	$coll_lineage				=	$homolka_hash->{$homolka_lineage} if(exists $homolka_hash->{$homolka_lineage});
	return($coll_lineage);
}



sub translate_coll2homolka {
	my $coll_lineage			=	shift;
	my $homolka_lineage                     =       'unkown';
	my $coll_hash				=	{};
	$coll_hash->{'1'}			=	'EAI';				# original	East-African-Indian
	$coll_hash->{'1.1'}			= 	'EAI';				# original	East-African-Indian
	$coll_hash->{'1.1.1'}			= 	'EAI';				# original	East-African-Indian
	$coll_hash->{'1.1.1.1'}			= 	'EAI'; 				# original	East-African-Indian
	$coll_hash->{'1.1.2'}			= 	'EAI';				# original	East-African-Indian
	$coll_hash->{'1.1.3'}			= 	'EAI';				# original	East-African-Indian
	$coll_hash->{'1.2.1'}			= 	'EAI Manila';			# original	East-African-Indian
	$coll_hash->{'1.2.2'}			= 	'EAI';				# original	East-African-Indian
	$coll_hash->{'2'}			= 	'East-Asian';
	$coll_hash->{'2.1'}			= 	'East-Asian non-Beijing';	# original East-Asian (non-Beijing)
	$coll_hash->{'2.2'}			= 	'Beijing';
	$coll_hash->{'2.2.1'}			= 	'Beijing';
	$coll_hash->{'2.2.1.1'}			= 	'Beijing'; 
	$coll_hash->{'2.2.1.2'}			= 	'Beijing';
	$coll_hash->{'2.2.2'}			= 	'Beijing';	
	$coll_hash->{'3'}			= 	'Delhi-CAS';			# original Indo-Oceanic
	$coll_hash->{'3.1.1'}			= 	'Delhi-CAS';			# original Indo-Oceanic
	$coll_hash->{'3.1.2'}			= 	'Delhi-CAS';			# original Indo-Oceanic
	$coll_hash->{'3.1.2.1'}			= 	'Delhi-CAS';			# original Indo-Oceanic
	$coll_hash->{'3.1.2.2'}			= 	'Delhi-CAS';			# original Indo-Oceanic
	$coll_hash->{'4'} 			= 	'Euro-American';
	$coll_hash->{'4.1'}			= 	'Euro-American';
	$coll_hash->{'4.1.1'}			= 	'X-type';			# original Euro-American (X-type)
	$coll_hash->{'4.1.1.1'}			= 	'X-type';			# original Euro-American (X-type)
	$coll_hash->{'4.1.1.2'}			= 	'X-type';			# original Euro-American (X-type)
	$coll_hash->{'4.1.1.3'}			= 	'X-type';			# original Euro-American (X-type)
	$coll_hash->{'4.1.2'}			= 	'Euro-American';
	$coll_hash->{'4.1.2.1'}			= 	'Haarlem';			# original Euro-American (Haarlem)
	$coll_hash->{'4.2'}			= 	'Euro-American';
	$coll_hash->{'4.2.1'}			= 	'Ural';				# original Euro-American (Ural)
	$coll_hash->{'4.2.2'}			= 	'Euro-American';
	$coll_hash->{'4.2.2.1'}			= 	'TUR';				# original Euro-American (TUR)
	$coll_hash->{'4.3'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.1'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.2'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.2.1'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.3'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.4'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.4.1'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.4.2'}			= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.3.4.2.1'}		= 	'LAM';				# original Euro-American (LAM)
	$coll_hash->{'4.4'}			= 	'Euro-American';
	$coll_hash->{'4.4.1'}			= 	'Euro-American';
	$coll_hash->{'4.4.1.1'}			= 	'S-type';			# original Euro-American (S-type)
	$coll_hash->{'4.4.1.2'}			= 	'Euro-American';
	$coll_hash->{'4.4.2'}			= 	'Euro-American';
	$coll_hash->{'4.5'}			= 	'Euro-American';
	$coll_hash->{'4.6'}			= 	'Euro-American';
	$coll_hash->{'4.6.1'}			= 	'Uganda';			# original Euro-American (Uganda)
	$coll_hash->{'4.6.1.1'}			= 	'Uganda';			# original Euro-American (Uganda)
	$coll_hash->{'4.6.1.2'}			= 	'Uganda';			# original Euro-American (Uganda)
	$coll_hash->{'4.6.2'}			= 	'Euro-American';
	$coll_hash->{'4.6.2.1'}			= 	'Euro-American';
	$coll_hash->{'4.6.2.2'}			= 	'Cameroon';			# original Euro-American (Cameroon)
	$coll_hash->{'4.7'}			= 	'mainly T';			# original Euro-American (mainly T)
	$coll_hash->{'4.8'}			= 	'mainly T';			# original Euro-American (mainly T)
	$coll_hash->{'4.9'}			= 	'H37Rv-like';			# original Euro-American (H37Rv-like)
	$coll_hash->{'5'} 			= 	'West-Africa 1';
	$coll_hash->{'6'} 			= 	'West-Africa 2';
	$coll_hash->{'7'} 			= 	'Ethiopia';
	$coll_hash->{'BOV'}			= 	'M. bovis';
	$homolka_lineage 			=	$coll_hash->{$coll_lineage} if (exists $coll_hash->{$coll_lineage});
	return($homolka_lineage);
}



#############################################################
###                                                       ###
###                     Specificator                      ###
###                                                       ###
#############################################################



sub specificator_beijing_easy {
	my $IDpositions				=	shift;;
	my $species				=	'unknown';
	my $lineage				=	'unknown';
	if($IDpositions->{782634}	eq "G") { 	$lineage = 'Asian/Africa';			}	#Beijing Subgroup Asian/Africa 1			
	if($IDpositions->{3061703}	eq "A") {	$lineage = 'Asian/Africa';			}	#Beijing Subgroup Asian/Africa 2			
	if($IDpositions->{2428517}	eq "G") {	$lineage = 'Pacific RD150';			}	#Beijing Subgroup Pacific RD150	
	if($IDpositions->{1584762}	eq "A") {	$lineage = 'Europe/Russian W148 Outbreak';	}	#Beijing Subgroup Europe/Russian W148 Outbreak	
	if($IDpositions->{4238675}	eq "T") {	$lineage = 'Central Asia';			}	#Beijing Subgroup Central Asia
	if($IDpositions->{46324}	eq "C") {	$lineage = 'Central Asia outbreak';		}	#Beijing Subgroup Central Asia outbreak
	if($IDpositions->{15890}	eq "A") {	$lineage = 'Ancestral 1'; 			}	#Beijing Subgroup Ancestral 1
	if($IDpositions->{4153340}	eq "A") {	$lineage = 'Ancestral 2';			}	#Beijing Subgroup Ancestral 2
	if($IDpositions->{3373397}	eq "A") {	$lineage = 'Ancestral 3';			}	#Beijing Subgroup Ancestral 3
	return($species,$lineage);
}	



sub specificator_coll_easy {
	my $IDpositions			=		shift;
	my $species			=		'unknown';
	my $lineage			=		'unknown';
	if($IDpositions->{931123}	eq "T")		{ $lineage = '4'		unless(length($lineage) > 1);	}		# wt allel L4 Euro-American
	if($IDpositions->{62657} 	eq "A")  	{ $lineage = '4.1' 		unless(length($lineage) > 3);	}		# L4.1 Euro-American
	if($IDpositions->{514245} 	eq "T") 	{ $lineage = '4.1.1' 		unless(length($lineage) > 5);	}		# L4.1.1 Euro-American (X-type)
	if($IDpositions->{1850119} 	eq "T") 	{ $lineage = '4.1.1.1' 		unless(length($lineage) > 7);	}		# L4.1.1.1 Euro-American (X-type)
	if($IDpositions->{541048} 	eq "G") 	{ $lineage = '4.1.1.2' 		unless(length($lineage) > 7);	}		# L4.1.1.2 Euro-American (X-type)
	if($IDpositions->{4229087} 	eq "T") 	{ $lineage = '4.1.1.3' 		unless(length($lineage) > 7);	}		# L4.1.1.3 Euro-American (X-type)
	if($IDpositions->{891756} 	eq "G") 	{ $lineage = '4.1.2' 		unless(length($lineage) > 5);	}		# L4.1.2 Euro-American
	if($IDpositions->{107794} 	eq "T") 	{ $lineage = '4.1.2.1' 		unless(length($lineage) > 7);	}		# L4.1.2.1 Euro-American (Haarlem)
	if($IDpositions->{2411730} 	eq "C")  	{ $lineage = '4.2' 		unless(length($lineage) > 3);	}		# L4.2  Euro-American
	if($IDpositions->{783601} 	eq "C") 	{ $lineage = '4.2.1' 		unless(length($lineage) > 5);	}		# L4.2.1 Euro-American (Ural)
	if($IDpositions->{1487796} 	eq "A") 	{ $lineage = '4.2.2' 		unless(length($lineage) > 5);	}		# L4.2.2  Euro-American
	if($IDpositions->{1455780} 	eq "C") 	{ $lineage = '4.2.2.1' 		unless(length($lineage) > 7);	}		# L4.2.2.1 Euro-American (TUR)
	if($IDpositions->{764995} 	eq "G")   	{ $lineage = '4.3' 		unless(length($lineage) > 3);	}		# L4.3 Euro-American (LAM)
	if($IDpositions->{615614} 	eq "A") 	{ $lineage = '4.3.1' 		unless(length($lineage) > 5);	}		# L4.3.1 Euro-American (LAM)
	if($IDpositions->{4316114} 	eq "A") 	{ $lineage = '4.3.2' 		unless(length($lineage) > 5);	}		# L4.3.2 Euro-American (LAM)
	if($IDpositions->{3388166} 	eq "G") 	{ $lineage = '4.3.2.1' 		unless(length($lineage) > 7);	}		# L4.3.2.1 Euro-American (LAM)
	if($IDpositions->{403364} 	eq "A") 	{ $lineage = '4.3.3' 		unless(length($lineage) > 5);	}		# L4.3.3 Euro-American (LAM)
	if($IDpositions->{3977226} 	eq "A") 	{ $lineage = '4.3.4' 		unless(length($lineage) > 5);	}		# L4.3.4 Euro-American (LAM)
	if($IDpositions->{4398141} 	eq "A") 	{ $lineage = '4.3.4.1' 		unless(length($lineage) > 7);	}		# L4.3.4.1 Euro-American (LAM)
	if($IDpositions->{1132368} 	eq "T") 	{ $lineage = '4.3.4.2' 		unless(length($lineage) > 7);	}		# L4.3.4.2 Euro-American (LAM)
	if($IDpositions->{1502120} 	eq "A") 	{ $lineage = '4.3.4.2.1' 	unless(length($lineage) > 9);	}		# L4.3.4.2.1 Euro-American (LAM)
	if($IDpositions->{4307886} 	eq "A")  	{ $lineage = '4.4' 		unless(length($lineage) > 3);	}		# L4.4  Euro-American
	if($IDpositions->{4151558} 	eq "A") 	{ $lineage = '4.4.1' 		unless(length($lineage) > 5);	}		# L4.4.1  Euro-American
	if($IDpositions->{355181} 	eq "A") 	{ $lineage = '4.4.1.1' 		unless(length($lineage) > 7);	}		# L4.4.1.1 Euro-American (S-type)
	if($IDpositions->{2694560} 	eq "C") 	{ $lineage = '4.4.1.2' 		unless(length($lineage) > 7);	}		# L4.4.1.2 Euro-American
	if($IDpositions->{4246508} 	eq "A") 	{ $lineage = '4.4.2' 		unless(length($lineage) > 5);	}		# L4.4.2  Euro-American
	if($IDpositions->{1719757} 	eq "T")  	{ $lineage = '4.5' 		unless(length($lineage) > 3);	}		# L4.5  Euro-American
	if($IDpositions->{3466426} 	eq "A")  	{ $lineage = '4.6' 		unless(length($lineage) > 3);	}		# L4.6  Euro-American
	if($IDpositions->{4260268} 	eq "C") 	{ $lineage = '4.6.1' 		unless(length($lineage) > 5);	}		# L4.6.1 Euro-American (Uganda)
	if($IDpositions->{874787} 	eq "A") 	{ $lineage = '4.6.1.1' 		unless(length($lineage) > 7);	}		# L4.6.1.1  Euro-American (Uganda)
	if($IDpositions->{1501468} 	eq "C") 	{ $lineage = '4.6.1.2' 		unless(length($lineage) > 7);	}		# L4.6.1.2  Euro-American (Uganda)
	if($IDpositions->{4125058} 	eq "C") 	{ $lineage = '4.6.2' 		unless(length($lineage) > 5);	}		# L4.6.2  Euro-American
	if($IDpositions->{3570528} 	eq "G") 	{ $lineage = '4.6.2.1' 		unless(length($lineage) > 7);	}		# L4.6.2.1  Euro-American
	if($IDpositions->{2875883} 	eq "T") 	{ $lineage = '4.6.2.2' 		unless(length($lineage) > 7);	}		# L4.6.2.2 Euro-American (Cameroon)
	if($IDpositions->{4249732} 	eq "G")  	{ $lineage = '4.7' 		unless(length($lineage) > 3);	}		# L4.7 Euro-American (mainly T)
	if($IDpositions->{3836739} 	eq "A")  	{ $lineage = '4.8' 		unless(length($lineage) > 3);	}		# L4.8 Euro-American (mainly T)
	if($IDpositions->{1759252} 	eq "G")  	{ $lineage = '4.9' 		unless(length($lineage) > 3);	}		# wt allel L4.9 Euro-American (H37Rv-like)
	if($IDpositions->{615938} 	eq "A")  	{ $lineage = '1'		unless(length($lineage) > 1);	}		# L1 East-African-Indian
	if($IDpositions->{4404247} 	eq "A")  	{ $lineage = '1.1'		unless(length($lineage) > 3);	}		# L1.1 East-African-Indian
	if($IDpositions->{3021283} 	eq "A")  	{ $lineage = '1.1.1'		unless(length($lineage) > 5);	}		# L1.1.1 East-African-Indian
	if($IDpositions->{3216553} 	eq "A")  	{ $lineage = '1.1.1.1'		unless(length($lineage) > 7);	}		# L1.1.1.1 East-African-Indian
	if($IDpositions->{2622402}	eq "A") 	{ $lineage = '1.1.2'		unless(length($lineage) > 5);	}		# L1.1.2 East-African-Indian
	if($IDpositions->{1491275} 	eq "A") 	{ $lineage = '1.1.3'		unless(length($lineage) > 5);	}		# L1.1.3 East-African-Indian
	if($IDpositions->{3479545} 	eq "A")  	{ $lineage = '1.2.1'		unless(length($lineage) > 5);	}		# L1.2.1 East-African-Indian
	if($IDpositions->{3470377} 	eq "T")  	{ $lineage = '1.2.2'		unless(length($lineage) > 5);	}		# L1.2.2 East-African-Indian
	if($IDpositions->{497491} 	eq "A") 	{ $lineage = '2'		unless(length($lineage) > 1);	}		# L2 East-Asian
	if($IDpositions->{1881090} 	eq "T")  	{ $lineage = '2.1'		unless(length($lineage) > 3);	}		# L2.1 East-Asian (non-Beijing)
	if($IDpositions->{2505085} 	eq "A")  	{ $lineage = '2.2'		unless(length($lineage) > 3);	}		# L2.2 East-Asian (Beijing)
	if($IDpositions->{797736} 	eq "T")  	{ $lineage = '2.2.1'		unless(length($lineage) > 5);	}		# L2.2.1 East-Asian (Beijing)
	if($IDpositions->{4248115} 	eq "T")  	{ $lineage = '2.2.1.1'		unless(length($lineage) > 7);	}		# L2.2.1.1 East-Asian (Beijing)
	if($IDpositions->{3836274} 	eq "A") 	{ $lineage = '2.2.1.2'		unless(length($lineage) > 7);	}		# L2.2.1.2 East-Asian (Beijing)
	if($IDpositions->{346693} 	eq "T")  	{ $lineage = '2.2.2'		unless(length($lineage) > 5);	}		# L2.2.2 East-Asian (Beijing)
	if($IDpositions->{3273107} 	eq "A") 	{ $lineage = '3'		unless(length($lineage) > 1);	}		# L3 Indo-Oceanic
	if($IDpositions->{1084911} 	eq "A") 	{ $lineage = '3.1.1'		unless(length($lineage) > 5);	}		# L3.1.1 Indo-Oceanic
	if($IDpositions->{3722702} 	eq "C") 	{ $lineage = '3.1.2'		unless(length($lineage) > 5);	}		# L3.1.2 Indo-Oceanic
	if($IDpositions->{1237818} 	eq "G")  	{ $lineage = '3.1.2.1'		unless(length($lineage) > 7);	}		# L3.1.2.1 Indo-Oceanic
	if($IDpositions->{2874344} 	eq "A")  	{ $lineage = '3.1.2.2'		unless(length($lineage) > 7);	}		# L3.1.2.2 Indo-Oceanic
	if($IDpositions->{1799921} 	eq "A") 	{ $lineage = '5'		unless(length($lineage) > 1);	}		# L5 West-Africa 1
	if($IDpositions->{1816587} 	eq "G") 	{ $lineage = '6'		unless(length($lineage) > 1);	}		# L6 West-Africa 2
	if($IDpositions->{1137518} 	eq "A") 	{ $lineage = '7'		unless(length($lineage) > 1);	}		# L7 Lineage 7
	if($IDpositions->{2831482} 	eq "G") 	{ $lineage = 'BOV'		unless(length($lineage) > 1);	}		# LBOV M. bovis
	return($species,$lineage);
}



sub specificator_coll_branch {
	my $IDpositions		=	shift;
	my $species		=	'unknown';
	my $lineage 		=	'unknown';
	if($IDpositions->{931123} eq "T")  {								# wt allel L4 Euro-American
		if($IDpositions->{1759252} eq "G")  {							# wt allel L4.9 Euro-American (H37Rv-like)
			$lineage = '4.9';
		}
		elsif($IDpositions->{1759252} eq "T") {							# L4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8
			if($IDpositions->{62657} eq "A")  {						# L4.1 Euro-American
				if($IDpositions->{514245} eq "T") {					# L4.1.1 Euro-American (X-type)
					if($IDpositions->{1850119} eq "T") {				# L4.1.1.1 Euro-American (X-type)
						$lineage = '4.1.1.1';
					}
					elsif($IDpositions->{541048} eq "G") {				# L4.1.1.2 Euro-American (X-type)
						$lineage = '4.1.1.2';
					}
					elsif($IDpositions->{4229087} eq "T") {				# L4.1.1.3 Euro-American (X-type)
						$lineage = '4.1.1.3';
					}
					else {
						$lineage = '4.1.1';
					}
				}
				elsif($IDpositions->{891756} eq "G") {					# L4.1.2 Euro-American
					if($IDpositions->{107794} eq "T") {				# L4.1.2.1 Euro-American (Haarlem)
						$lineage = '4.1.2.1';
					}
					else {
						$lineage = '4.1.2';
					}
				}
				else {
					$lineage = '4.1';
				}

			}
			elsif($IDpositions->{2411730} eq "C")  {					# L4.2  Euro-American
				if($IDpositions->{783601} eq "C") {					# L4.2.1 Euro-American (Ural)
					$lineage = '4.2.1';
				}
				elsif($IDpositions->{1487796} eq "A") {					# L4.2.2  Euro-American
					if($IDpositions->{1455780} eq "C") {				# L4.2.2.1 Euro-American (TUR)
						$lineage = '4.2.2.1';
					}
					else {
						$lineage = '4.2.2';
					}
				}
				else {
					$lineage = '4.2';
				}
			}
			elsif($IDpositions->{764995} eq "G")   {					# L4.3 Euro-American (LAM)
				if($IDpositions->{615614} eq "A") {					# L4.3.1 Euro-American (LAM)
					$lineage = '4.3.1';
				}
				elsif($IDpositions->{4316114} eq "A") {					# L4.3.2 Euro-American (LAM)
					if($IDpositions->{3388166} eq "G") {				# L4.3.2.1 Euro-American (LAM)
						$lineage = '4.3.2.1';
					}
					else {
						$lineage = '4.3.2';
					}
				}
				elsif($IDpositions->{403364} eq "A") {					# L4.3.3 Euro-American (LAM)
					$lineage = '4.3.3';
				}
				elsif($IDpositions->{3977226} eq "A") {					# L4.3.4 Euro-American (LAM)
					if ($IDpositions->{4398141} eq "A") {				# L4.3.4.1 Euro-American (LAM)
						$lineage = '4.3.4.1';
					}
					elsif($IDpositions->{1132368} eq "T") {				# L4.3.4.2 Euro-American (LAM)
						if($IDpositions->{1502120} eq "A") {			# L4.3.4.2.1 Euro-American (LAM)
							$lineage = '4.3.4.2.1';
						}
						else {
							$lineage = '4.3.4.2';
						}
					}
					else {
						$lineage = '4.3.4';
					}
				}
				else {
					$lineage = '4.3';
				}
			}
			elsif($IDpositions->{4307886} eq "A")  {					# L4.4  Euro-American
				if($IDpositions->{4151558} eq "A") {					# L4.4.1  Euro-American
					if($IDpositions->{355181} eq "A") {				# L4.4.1.1 Euro-American (S-type)
						$lineage = '4.4.1.1';

					}
					elsif($IDpositions->{2694560} eq "C") {				# L4.4.1.2 Euro-American
						$lineage = '4.4.1.2';
					}
					else {
						$lineage = '4.4.1';
					}
				}
				elsif($IDpositions->{4246508} eq "A") {					# L4.4.2  Euro-American
					$lineage = '4.4.2';
				}
				else {
					$lineage = '4.4';
				}
			}
			elsif($IDpositions->{1719757} eq "T")  {					# L4.5  Euro-American
				$lineage = '4.5';
			}
			elsif($IDpositions->{3466426} eq "A")  {					# L4.6  Euro-American
				if ($IDpositions->{4260268} eq "C") {					# L4.6.1 Euro-American (Uganda)
					if ($IDpositions->{874787} eq "A") {				# L4.6.1.1  Euro-American (Uganda)
						$lineage = '4.6.1.1';
					}
					elsif ($IDpositions->{1501468} eq "C") {			# L4.6.1.2  Euro-American (Uganda)
						$lineage = '4.6.1.2';
					}
					else {
						$lineage = '4.6.1';
					}
				}
				elsif($IDpositions->{4125058} eq "C") {					# L4.6.2  Euro-American
					if($IDpositions->{3570528} eq "G") {				# L4.6.2.1  Euro-American
						$lineage = '4.6.2.1';
					}
					elsif($IDpositions->{2875883} eq "T") {				# L4.6.2.2 Euro-American (Cameroon)
						$lineage = '4.6.2.2';
					}
					else {
						$lineage = '4.6.2';
					}
				}
				else {
					$lineage = '4.6';
				}
			}
			elsif($IDpositions->{4249732} eq "G")  {					# L4.7 Euro-American (mainly T)
				$lineage = '4.7';
			}
			elsif($IDpositions->{3836739} eq "A")  {					# L4.8 Euro-American (mainly T)
				$lineage = '4.8';
			}
		}
		else {
			$lineage = '4';
		}
	}
	else {												# L1,2,3,5,6,7,BOV,BOV_AFRI	
		if($IDpositions->{615938} eq "A")  { 							# L1 Indo-Oceanic
			if($IDpositions->{4404247} eq "A")  {						# L1.1 Indo-Oceanic
				if($IDpositions->{3021283} eq "A")  {					# L1.1.1 Indo-Oceanic
					if($IDpositions->{3216553} eq "A")  {				# L1.1.1.1 Indo-Oceanic
						$lineage = '1.1.1.1'; 
					} 
					else {
						$lineage = '1.1.1';
					}
				}
				elsif($IDpositions->{2622402} eq "A") {					# L1.1.2 Indo-Oceanic
					$lineage = '1.1.2';
				}
				elsif($IDpositions->{1491275} eq "A") {					# L1.1.3 Indo-Oceanic
					$lineage = '1.1.3';
				}
				else {
				$lineage = '1.1';
				}
			}
			elsif($IDpositions->{3479545} eq "A")  {					# L1.2.1 Indo-Oceanic
				$lineage = '1.2.1';
				}
			elsif($IDpositions->{3470377} eq "T")  {					# L1.2.2 Indo-Oceanic
				$lineage = '1.2.2';
				}
			else {
				$lineage = '1';
			}
		}
			
		elsif($IDpositions->{497491} eq "A") {							# L2 East-Asian
			if($IDpositions->{1881090} eq "T")  {
				$lineage = '2.1';							# L2.1 East-Asian (non-Beijing)
			}
			elsif($IDpositions->{2505085} eq "A")  {					# L2.2 East-Asian (Beijing)
				if($IDpositions->{797736} eq "T")  {					# L2.2.1 East-Asian (Beijing)
					if($IDpositions->{4248115} eq "T")  {				# L2.2.1.1 East-Asian (Beijing)
						$lineage = '2.2.1.1'; 
					}
 					elsif($IDpositions->{3836274} eq "A") {				# L2.2.1.2 East-Asian (Beijing)
						$lineage = '2.2.1.2';
					}
					else {
						$lineage = '2.2.1';
					}
				}
				elsif($IDpositions->{346693} eq "T")  {					# L2.2.2 East-Asian (Beijing)
					$lineage = '2.2.2';
				}
				else {
					$lineage = '2.2';
				}
			}
			else {
				$lineage = '2';
			}
		}
		elsif($IDpositions->{3273107} eq "A") {							# L3 East-African-Indian
			if($IDpositions->{1084911} eq "A") {						# L3.1.1 East-African-Indian
				$lineage = '3.1.1';
			}
			elsif($IDpositions->{3722702} eq "C") {						# L3.1.2 East-African-Indian
				if ($IDpositions->{1237818} eq "G")  {					# L3.1.2.1 East-African-Indian
					$lineage = '3.1.2.1';
				}
				elsif ($IDpositions->{2874344} eq "A")  {				# L3.1.2.2 East-African-Indian
					$lineage = '3.1.2.2';
				}
				else {
					$lineage = '3.1.2';
				}	
			}
			else {
				$lineage = '3';
			}
		}
		elsif ($IDpositions->{1799921} eq "A") {						#L5 West-Africa 1
			$lineage = '5';
		}
		elsif ($IDpositions->{1816587} eq "G") {						#L6 West-Africa 2
			$lineage = '6';
		}
		elsif ($IDpositions->{1137518} eq "A") {						#L7 Lineage 7
			$lineage = '7';
		}
		elsif ($IDpositions->{2831482} eq "G") {						#LBOV M. bovis
			$lineage = 'BOV';
		}
	}
	return($species,$lineage);
}



sub specificator_homolka {
	my $IDpositions		= 	shift;
	my $species		= 	'unknown';
	my $lineage		=	'unknown';
	if($IDpositions->{648856} eq 'C') {									# not Euro American
		if($IDpositions->{648756} eq 'T') {								# M. africanum 1a/1b
			if($IDpositions->{2053726} eq 'T') {							# M. africanum 1b
				$species = 'M. africanum'; 
				$lineage = 'West African 1b'; 
			} 
			else {											# M. africanum 1a ???
				$species = 'M. africanum';
				$lineage = 'West African 1a';							# spoke with susanne, if 1b SNP is not there, it's 1a
			}
		}
		elsif($IDpositions->{649345} eq 'T') {								# EAI
			if($IDpositions->{1128814} eq 'A') {							# EAI Manila
				$species = 'M. tuberculosis';
				$lineage = 'EAI Manila'; 
			}
			else {											# EAI
				$species = 'M. tuberculosis';
				$lineage = 'EAI'; 
			}
		}
		elsif($IDpositions->{649585} eq 'A') {								# M. caprae
			$species = 'M. caprae';
		}
		elsif($IDpositions->{649601} eq 'T') {								# M. canettii
			$species = 'M. canettii'; 
		}
		else {												# Rv0557 is WT
			if($IDpositions->{157129} eq 'T') {     						# Delhi/CAS
				$species = 'M. tuberculosis';
				$lineage = 'Delhi-CAS'; 
			}
			elsif($IDpositions->{1128825} eq 'T') {     							# M. microti or M. pinnipedii
				$species = 'M. microti or M. pinnipedii'; 
				if($IDpositions->{1473079} eq 'A') {							# M. microti (not SH)
					$species = 'M. microti'; 
				}
				elsif(($IDpositions->{1674520} eq 'T') || ($IDpositions->{1473094} eq 'C')) {		# M. pinnipedii (not SH)
					$species = 'M. pinnipedii'; 
				}
			}
			elsif($IDpositions->{1129160} eq 'T') {     							# M. bovis
				$species = 'M. bovis';
			}
			elsif($IDpositions->{2053762} eq 'T') {     							# M. africanum 2
				$species = 'M. africanum';
				$lineage = 'West African 2'; 
			}
			elsif($IDpositions->{2955957} eq 'C') {     							# Beijing
     				$species = 'M. tuberculosis';
				$lineage = 'Beijing';
			} 
		} 
	}
	else {														# Euro American
		if($IDpositions->{648990} eq 'C') {									# Haarlem
			$species = 'M. tuberculosis';
			$lineage = 'Haarlem'; 
		}
		elsif($IDpositions->{648992} eq 'G') {									# S-type
			$species = 'M. tuberculosis';
			$lineage = 'S-type'; 
		}
		elsif($IDpositions->{649067} eq 'G') {									# Cameroon
			$species = 'M. tuberculosis';
			$lineage = 'Cameroon'; 
		}
		else {													# Rv0557 is WT
			if($IDpositions->{157292} eq 'T') {   	   							# LAM
				$species = 'M. tuberculosis';
				$lineage = 'LAM'; 
			}
			elsif($IDpositions->{1129165} eq 'A') {	     							# TUR
				$species = 'M. tuberculosis';
				$lineage = 'TUR';
			}
			elsif($IDpositions->{2053454} eq 'A') {	     							# Ural
				$species = 'M. tuberculosis';
				$lineage = 'Ural'; 
			}
			elsif($IDpositions->{2956731} eq 'T') {	     							# Ghana
				$species = 'M. tuberculosis';
				$lineage = 'Ghana';
			}
			elsif($IDpositions->{7539} eq 'G') {		     						# Uganda (not SH)
				$species = 'M. tuberculosis';
				$lineage = 'Uganda'; 
			}
		}
	}
	if($species eq 'unknown') {
		if($IDpositions->{1816848} eq 'T') {
			$species = 'M. africanum or M. caprae or M. microti or M. pinipedii'; 
		}		
		elsif($IDpositions->{2955233} eq 'T') {
			$species = 'M. africanum or M.bovis or M. canetii or M. caprae or M. microti or M. pinipedii or M. tuberculosis EAI';
		}
		elsif(!($IDpositions->{2955233} eq 'T')) {
			$species = 'M. tuberculosis';
			$lineage = 'Clade 1';
		}
	}
	return($species,$lineage);
}


1;
