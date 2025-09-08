#!/usr/bin/env perl

package TBresi;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

##needed for resi --> add to requirements!
use Statistics::R;
use Getopt::Long;
use File::Basename;
use File::Slurp;

$VERSION    =  1.0.0;
@ISA        =  qw(Exporter);
@EXPORT     =  qw(tbresi
                  tbresisummary
                  tbcombinedresi);


sub tbresi {
   # get parameter and input from front-end.
   my $logprint                  =  shift;
   my $VAR_dir                   =  shift;
   my $CALL_OUT                  =  shift;
   my $RESI_OUT                  =  shift;
   my $date_string               =  shift;
   my $resi_list_master          =  shift;
   my $resi_list_date            =  shift;
   my @truecodon_files           =     @_;
   my $minion                    =     "";
   my $res_hash                  =     {};
   my $change_hash               =     {};
   
   print $logprint ("<INFO>\t",timer(),"\tNo resistance file $resi_list_master. Will skip resistance annotation.\n") unless(-f $resi_list_master);
   #if($resi_list_master            eq   ''||"NONE"   ) { die "\n[ERROR]\t",timer(),"\tNeed to provide a resistance file. Use --help for usage information\n";}

   ###aminoacid changes###
   my %codon = ('Ser'=>'S',
			'Phe'=>'F',
			'Leu'=>'L',
			'Tyr'=>'Y',
			'_'=>'_',
			'Cys'=>'C',
			'Trp'=>'W',
			'Leu'=>'L',
			'Pro'=>'P',
			'His'=>'H',
			'Gln'=>'Q',
			'Arg'=>'R',
			'Ile'=>'I',
			'Met'=>'M',
			'Thr'=>'T',
			'Asn'=>'N',
			'Lys'=>'K',
			'Val'=>'V',
			'Ala'=>'A',
			'Asp'=>'D',
			'Glu'=>'E',
			'Gly'=>'G'
			);
   print $logprint "<INFO>\t",timer(),"\tStart parsing $resi_list_master...\n";
   if ($resi_list_master ne ""||"NONE"){
   open(Fin, "<", $resi_list_master) or die "\n[ERROR]\t",timer(),"\tUnable to open $resi_list_master: TBresi.pm line: ", __LINE__ , " \n";
   <Fin>;
   while (my $line = <Fin>){
   $line          =~  s/\015?\012?$//;
   my @line = split("\t",$line);
   my $pos = $line[0];
   ($pos = $pos) =~s/[i|\+]\d+$//;
   my $ref_allel = $line[1];
   my $mut_allel = $line[2];
   #if ($mut_allel =~ /^INS_/){print "$mut_allel\n";}
   $mut_allel =~s/^INS_//;
   #if ($mut_allel =~ /A|G|T|C/){print "$mut_allel\n";}
   my $type = $line[3];
   my $affects = $line[4];
   my $gene = $line[5];
   my $orientation = $line[6];
   my $gene_name = $line[7];
   if ($gene_name eq "gidB"){$gene_name="gid"};
   my $change = $line[8];
   my $predict_phyl = $line[9];
   my $predict = $line[10];
   my $comment = $line[13];
   my $mic = $line[14];
   my $whotier = $line[19];

	if ($pos ne "NA"){
		$change_hash->{$pos}->{$mut_allel} =   {change       =>   $change,
												predict      =>   $predict,
												predict_phyl =>   $predict_phyl,
												affects      =>   $affects,
												mut_allel    =>   $mut_allel,
												gene_name    =>   $gene_name,
												gene         =>   $gene,
												comment      =>   $comment};
	}
	$res_hash->{$gene}->{$change} = {predict => $predict,
									predict_phyl => $predict_phyl,
									gene_name => $gene_name,
									affects => $affects,
									comment => $comment};
	}
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $resi_list_master!\n";
   close (Fin);
   }

   foreach my $variant_file (@truecodon_files) {
      print $logprint ("<INFO>\t",timer(),"\tNo variant file $variant_file. Will skip resistance annotation.\n") unless(-f "$CALL_OUT/$variant_file");
      print $logprint ("<INFO>\t",timer(),"\tStart calling resistance for $variant_file...\n");
		if ($variant_file =~ /.*gatk_position_true-codon-variants.*/){
			(my $out_file = basename($variant_file)) =~ s/\..*$//g;
			my $output_mode = $1 if ($variant_file) =~/^.*outmode(\d\d\d)/;

			open(OUT,">","$RESI_OUT/${out_file}.gatk_position_true-codon-variants_outmode${output_mode}_${resi_list_date}_resi.tab") or die "\n<ERROR>\t",timer(),"\tUnable to create $${out_file}.gatk_position_true-codon-variants_outmode${output_mode}_${resi_list_date}_resi.tab: TBresi.pm line: ", __LINE__ , " \n";


			open (IN, "<", "$CALL_OUT/$variant_file") or die "\n[ERROR]\t",timer(),"\tUnable to open $variant_file: TBresi.pm line: ", __LINE__ , " \n";
			my $header = <IN>;
			$header =~ s/\015?\012?$//;
			$header =~ s/\tResistanceSNP\tPhyloSNP//;
			$header .= "\tResistance\tBenign\tMutation\tComment";
			print OUT "$header\n";

			while (my $line = <IN>) {
			chomp($line);
			next unless $line;
			next if $line =~ /^[A-Za-z]/;
			$line =~ s/\015?\012?$//;
			my @line = split("\t", $line);
			my $pos			= $line[0];		$pos 			= 0 unless($pos);
			my $pos_insertion = $pos;
			$pos_insertion =~s/[i|\+]\d+$//;
			my $ref			= $line[1];		$ref 			= 0 unless($ref);
			my $type		= $line[2];		$type 			= 0 unless($type);
			next if ($minion and $type eq "Del");
			my $allel		= $line[3];		$allel 			= 0 unless($allel);
			my $cov_forward = $line[4];		$cov_forward 	= 0 unless($cov_forward);
			my $cov_reverse = $line[5];		$cov_reverse 	= 0 unless($cov_reverse);
			my $qual_20		= $line[6];		$qual_20 		= 0 unless($qual_20);
			my $freq1		= $line[7];		$freq1 			= 0 unless($freq1);
			my $coverage	= $line[8];		$coverage 		= 0 unless($coverage);
			my $subs		= $line[9];		$subs 			= 0 unless($subs);
			my $subst		= $subs;
			$subst=~s/\s.*$//;
			#print "$subst\n";
			if (($subst ne "-") and $subst !~ /^\s*$/){
				my $old_aa;
				my $new_aa;
					if (substr($subst, 0, 3) !~ /\d+/) {
						$old_aa = substr( $subst, 0, 3 );
						# print "$old_aa\n";
						$old_aa = $codon{$old_aa};
						#print "$old_aa\n";
						}
					else{$old_aa = substr($subst, 0, 1)}

					if (substr($subst, -3) !~ /\d+/){
						$new_aa = substr( $subst, -3);
						$new_aa = $codon{$new_aa};
					}
					else{$new_aa = substr($subst, -1)}

				$subst =~ s/[^0-9]//g;
				$subst = $old_aa.$subst.$new_aa;
			}
		my $gene		= $line[10];	$gene 			= 0 unless($gene);
		my $gene_name	= $line[11];	$gene_name 		= 0 unless($gene_name);
		if ($gene_name eq "-"){$gene_name = $gene};
		if ($gene_name eq "gidB"){$gene_name = "gid"};
		my $annotation	= $line[12];	$annotation 	= 0 unless($annotation);
		#my $resistance	= $line[13];	$resistance 	= 0 unless($resistance);
		#my $phylo		= $line[14];	$phylo 			= 0 unless($phylo);
		my $region		= $line[15];	$region 		= 0 unless($region);
		my $warning		= $line[16];	$warning 		= "-" unless($warning);

		my $gene_position;
		my $wt_allel = $ref;
		my $subst_allel = $allel;
		#print "$subst_allel\n";
		my $better_res = "-";
		my $better_res_change = "-";
		my $benigninfo = "-";
		my $res_comment = "-";

#if ($gene_name eq "rrl"){print Dumper(\$res_hash->{$gene_name});}
			if (exists($res_hash->{$gene}->{$subst})){
			#if ($gene_name eq "rrl"){print Dumper(\$res_hash->{$gene_name});}
				if ($res_hash->{$gene}->{$subst}->{predict} ne "benign" and $res_hash->{$gene}->{$subst}->{predict} ne "likely benign"){
					$better_res = $res_hash->{$gene}->{$subst}->{predict};
					$gene_name = $res_hash->{$gene}->{$subst}->{gene_name};
					$better_res_change = $gene_name." ".$subst;
					$res_comment = $res_hash->{$gene}->{$subst}->{comment};
				}
				elsif ($res_hash->{$gene}->{$subst}->{predict} eq "benign" or $res_hash->{$gene}->{$subst}->{predict} eq "likely benign"){
					$benigninfo = $res_hash->{$gene}->{$subst}->{predict};
					$gene_name = $res_hash->{$gene}->{$subst}->{gene_name};
					$better_res_change = $gene_name." ".$subst;
					$res_comment = $res_hash->{$gene}->{$subst}->{comment};
				}
			}

			elsif (exists($change_hash->{$pos_insertion}->{$subst_allel})){
				if ($change_hash->{$pos_insertion}->{$subst_allel}->{predict} ne "benign" and $change_hash->{$pos_insertion}->{$subst_allel}->{predict} ne "likely benign" and ($subst_allel=~ /^((?!NA).)*[A|T|C|G]$/gm or $subst_allel eq "GAP")){
					$better_res = $change_hash->{$pos_insertion}->{$subst_allel}->{predict};
					$gene_name = $change_hash->{$pos_insertion}->{$subst_allel}->{gene_name};
					$better_res_change = $gene_name." ".$change_hash->{$pos_insertion}->{$subst_allel}->{change};
					$res_comment = $change_hash->{$pos_insertion}->{$subst_allel}->{comment};
				}
				elsif ($change_hash->{$pos_insertion}->{$subst_allel}->{predict} eq "benign" or $change_hash->{$pos_insertion}->{$subst_allel}->{predict} eq "likely benign" and ($subst_allel=~ /^((?!NA).)*[A|T|C|G]$/gm or $subst_allel eq "GAP")){
					#if ($pos_insertion == 1474959){print Dumper(\$change_hash->{$pos_insertion});print "$subst_allel\n";print "$pos_insertion\n";}
					$benigninfo = $change_hash->{$pos_insertion}->{$subst_allel}->{predict};
					$gene_name = $change_hash->{$pos_insertion}->{$subst_allel}->{gene_name};
					$better_res_change = $gene_name." ".$change_hash->{$pos_insertion}->{$subst_allel}->{change};
					$res_comment = $change_hash->{$pos_insertion}->{$subst_allel}->{comment};
				}
			}
			###general rules for dealing with Del, Ins and stop codons in certain genes
			
			###pncA reverse
			elsif ($gene_name eq "pncA" and $type eq "Del"){
			$gene_position = (2289241-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res = "PZA-R";
			}
			elsif ($gene_name eq "pncA" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (2289241-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res = "PZA-R";
				}
			}
			elsif ($gene_name eq "pncA" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="PZA-R";
				$better_res_change = $subst;
				}
				elsif($subst=~/(^\w).*(\w$)/ and $1 ne $2){
				$better_res ="PZA-R#";
				$better_res_change = $gene_name." ".$subst;
				$res_comment ="called by global rule, not confirmed";
				}
			}
			###katG reverse
			elsif ($gene_name eq "katG" and $type eq "Del"){
			#$res = 1;
			$gene_position = (2156111-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="INH-R";
			}
			elsif ($gene_name eq "katG" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (2156111-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="INH-R";
				}
			}
			elsif ($gene_name eq "katG" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="INH-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			###ethA reverse
			elsif ($gene_name eq "ethA" and $type eq "Del"){
			$gene_position = (4327473-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="ETH-R";
			}
			elsif ($gene_name eq "ethA" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (4327473-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="ETH-R";
				}
			}
			elsif ($gene_name eq "ethA" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="ETH-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			##gidB reverse - gene name = gid in variant file
			elsif (($gene_name eq "gid" or $gene_name eq "gidB") and $type eq "Del"){
			$gene_position = (4408202-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="SM-R";
			}
			elsif (($gene_name eq "gid" or $gene_name eq "gidB") and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (4408202-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="SM-R";
				}
			}
			elsif (($gene_name eq "gid" or $gene_name eq "gidB") and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="SM-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## rpoB forward
			elsif ($gene_name eq "rpoB" and $type eq "Del"){
			$gene_position = ($pos-759807)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="RIF-R";
			}
			elsif ($gene_name eq "rpoB" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-759807)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="RIF-R";
				}
			}
			elsif ($gene_name eq "rpoB" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="RIF-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## ald forward
			elsif ($gene_name eq "ald" and $type eq "Del"){
			$gene_position = ($pos-3086820)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="CS-R";
			}
			elsif ($gene_name eq "ald" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-3086820)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="CS-R";
				}
			}
			elsif ($gene_name eq "ald" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="CS-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## Rv0678 forward
			elsif ($gene eq "Rv0678" and $type eq "Del"){
			$gene_position = ($pos-778990)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="BDQ-R,CFZ-R";
			}
			elsif ($gene eq "Rv0678" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-778990)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="BDQ-R,CFZ-R";
				}
			}
			elsif ($gene eq "Rv0678" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="BDQ-R,CFZ-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			
			## ddn  Rv3547 forward
			elsif ($gene_name eq "ddn" and $type eq "Del"){
			$gene_position = ($pos-3986844)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="DLM-R";
			}
			elsif ($gene_name eq "ddn" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-3986844)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="DLM-R";
				}
			}
			elsif ($gene_name eq "ddn" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="DLM-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## tlyA Rv1694 forward
			elsif ($gene_name eq "tlyA" and $type eq "Del"){
			$gene_position = ($pos-1917940)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="CAP-R";
			}
			elsif ($gene_name eq "tlyA" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-1917940)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="CAP-R";
				}
			}
			elsif ($gene_name eq "tlyA" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="CAP-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}

   print OUT "$pos\t$ref\t$type\t$allel\t$cov_forward\t$cov_reverse\t$qual_20\t$freq1\t$coverage\t$subs\t$gene\t$gene_name\t$annotation\t$region\t$warning\t$better_res\t$benigninfo\t$better_res_change\t$res_comment\n";
   
			}
   close (IN);
   close (OUT);
   

		}
	print $logprint "<INFO>\t",timer(),"\tFinished calling resistance for $variant_file!\n";
	}
$res_hash           =  {};
$change_hash        =  {};
@truecodon_files    =  ();
}

sub tbresisummary {
    my $logprint                 =      shift;
    my $RESI_OUT                 =      shift;
	my $snp_vars                 =      shift;
	my $resi_list_date           =      shift;
	my $cutoff                   =      shift;
	my $qcutoff                  =      shift;
	my $fcutoff                  =      shift;
	my $pcutoff                  =      shift;
	my $mincutoff                =      shift;
	my @resi_files               =      @_;
	my $resi_file                =       "";
	my $file                     =       "";
	my $empty					 =		 "";	
	my @line_number				 ;
	my $line                     =       "";
	my @ID                       ;
	my $outputfile               =       "";
	my @names                   ;
	my $name                    =       "";
	my $dim                    =       "";
	my $antibiotic                    =       "";
	my $prediction                    =       "";
	my %mutations;


#opendir(RESIDIR,"$RESI_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $RESI_OUT: MTBseq.pl line: ", __LINE__ ," \n";
#@resi_files        =  grep { $_ =~ /^\w.*\.gatk_position_true-codon-variants_outmode[01]${snp_vars}1\_$resi_list_date\_resi.tab$/ && -f "$RESI_OUT/$_"   }  readdir(RESIDIR);
#closedir(RESIDIR);

  foreach my $resi_file (@resi_files) { #check whether the file is empty
	
	  $empty = "full";
	  my @line_number = read_file("$RESI_OUT/$resi_file");
	  my $lines = @line_number;
	  if ($lines == 1){
		  $empty = "empty";
	  }
	  next if ($empty eq "empty");
	  
	$resi_file=~/^(.+).tab/ or die "\n<ERROR>\t",timer(),"Strange file format: $resi_file: TBresi.pm line: ", __LINE__ , " \n";
	my $file=$1;
	open(Fout,">$RESI_OUT/${file}_summary.tab") or die "\n<ERROR>\t",timer(),"Cannot write output file: TBresi.pm line: ", __LINE__ , " \n";
	print Fout "SampleID\tLibID\tINH\tFreq_INH\tRMP\tFreq_RMP\tSM\tFreq_SM\tEMB\tFreq_EMB\tPZA\tFreq_PZA\tMFX\tFreq_MFX\tLFX\tFreq_LFX\tCFZ\tFreq_CFZ\tKAN\tFreq_KAN\tAMK\tFreq_AMK\tCPR\tFreq_CPR\tETH/PTH\tFreq_ETH/PTH\tLZD\tFreq_LZD\tBDQ\tFreq_BDQ\tCS\tFreq_CS\tPAS\tFreq_PAS\tDLM\tFreq_DLM\tPrediction\n";
    @ID=split("_",$file);
    
    open (Fin,"<$RESI_OUT/$resi_file") or die "\n<ERROR>\t",timer(),"Cannot open $file: TBresi.pm line: ", __LINE__ , " \n";
    print $logprint ("<INFO>\t",timer(),"\tStart summarizing resistance for $resi_file...\n");
	my $R = Statistics::R->new();
       $R->startR;
       $R->send('library("readr")');#R-Library for reading in tab-seperated files
        
       $R->run(qq'table<-read_delim(file="${RESI_OUT}/${resi_file}", "\t", escape_double =F, trim_ws=T);'); #Read in Table with resistance-calls 
       @names=("INH-R (INH)","RIF-R (RMP)","SM-R (SM)","EMB-R (EMB)","PZA-R (PZA)","MFX-R (MFX)","LFX-R (LFX)","CFZ-R (CFZ)","KAN-R (KAN)","AMI-R (AMK)","CAP-R (CPR)","ETH-R (ETH)","LZD-R (LZD)","BDQ-R (BDQ)","CS-R (CS)","PAS-R (PAS)","DLM-R (DLM)");
            foreach $a (@names){#loop over all antibiotics
    
			$R->run(qq'name<-"$a"
                cutoff<-as.numeric("$cutoff")
                qcutoff<-as.numeric("$qcutoff")
                fcutoff<-as.numeric("$fcutoff")
				pcutoff<-as.numeric("$pcutoff")
				mincutoff<-as.numeric("$mincutoff")
                short<-as.character(unlist(strsplit("$a"," ",fixed=TRUE))[1])
                anti<-as.character(unlist(strsplit("$a"," ",fixed=TRUE))[2])
                antibiotic<-gsub("[()]","", anti)');
			$name=$R->get('short');
			$antibiotic=$R->get('antibiotic');
            
			print "Search for mutations mediating resistance to $antibiotic\n";   
			$R->run('table1<-table[which(table$CovFor >=cutoff & table$CovRev >=cutoff & table$Qual20 >=qcutoff & table$Freq >=fcutoff & table$Freq <=25 & table$Cov >= mincutoff & (table$Qual20/(table$CovRev + table$CovFor)*100) >= pcutoff),]
				table2<-table[which(table$CovFor >=cutoff & table$CovRev >=cutoff & table$Qual20 >=qcutoff & table$Freq >=fcutoff & table$Freq >25 & table$Cov >= mincutoff),]
				table3<-rbind(table1, table2)
                idx<-table3[grep(short,table3$Resistance),]
                count<-dim(idx)[1]'); #subtable with all lines, that have an entry in the Resistance column for the specific drug with "-R" and fullfilling the thresholds
			$dim=$R->get('count');
        

            if($dim >= 1){#check if there are any resistance calls
            $R->run(q'mut<-1
                     freq<-1
                if(dim(idx)[1]>1){
                    for(i in 1:dim(idx)[1]){
                        
                            mut[i]<-idx$Mutation[i];
                            freq[i]<-round(idx$Freq[i],digits=2);
                    }
                    if(length(mut)==2){mut<-paste(mut[1],mut[2],sep="; ")
                                       freq<-paste(freq[1],freq[2],sep="; ")
                    }else if (length(mut)==3){mut<-paste(mut[1],mut[2],mut[3],sep="; ")
                                              freq<-paste(freq[1],freq[2], freq[3],sep="; ")
                    }else if (length(mut)==4){mut<-paste(mut[1],mut[2],mut[3], mut[4],sep="; ")
                                            freq<-paste(freq[1],freq[2], freq[3], freq[4],sep="; ")
                    }else if (length(mut)==5){mut<-paste(mut[1],mut[2],mut[3],mut[4], mut[5],sep="; ")
                                              freq<-paste(freq[1],freq[2], freq[3], freq[4],freq[5],sep="; ")
                    }else if (length(mut)==6){mut<-paste(mut[1],mut[2],mut[3],mut[4], mut[5], mut[6],sep="; ")
                                            freq<-paste(freq[1],freq[2], freq[3], freq[4],freq[5],freq[6],sep="; ")
                    }else {mut<-paste(mut[1],mut[2],mut[3],mut[4], mut[5], mut[6],"and others",sep="; ")
                    freq<-paste(freq[1],freq[2],freq[3],freq[4], freq[5], freq[6],"and others",sep="; ")}
                }else{
                        mut<-idx$Mutation
                        freq<-round(idx$Freq, digits=2)
                    }');
            
            $mutations{$antibiotic}=$R->get('mut');
            $mutations{$antibiotic.'_Freq'}=$R->get('freq');
           
            }

            else{$mutations{$antibiotic}="-";
                $mutations{$antibiotic.'_Freq'}="-";}

        
        }
    if ($mutations{INH} ne "-" && $mutations{RMP} ne "-" && ($mutations{MFX} ne "-" || $mutations{LFX} ne "-") && ($mutations{BDQ} ne "-" || $mutations{LZD} ne "-" )){$prediction = "XDR";
    }elsif ($mutations{INH} ne "-" && $mutations{RMP} ne "-" && (($mutations{MFX} ne "-" || $mutations{LFX} ne "-"))){$prediction = "preXDR";
    }elsif ($mutations{INH} ne "-" && $mutations{RMP} ne "-"){$prediction = "MDR";
    }elsif ($mutations{RMP} ne "-") {$prediction = "RR";
    }elsif ($mutations{INH} eq "-" && $mutations{RMP} eq "-" && $mutations{SM} eq "-" && $mutations{EMB} eq "-" && $mutations{PZA} eq "-" && $mutations{MFX} eq "-" && $mutations{LFX} eq "-" && $mutations{CFZ} eq "-" && $mutations{KAN} eq "-" && $mutations{AMK} eq "-" && $mutations{CPR} eq "-" && $mutations{ETH} eq "-" && $mutations{LZD} eq "-" && $mutations{BDQ} eq "-" && $mutations{CS} eq "-" && $mutations{PAS} eq "-" && $mutations{DLM} eq "-"){$prediction = "S";
    }else{$prediction = "nonMDR";}
    print Fout "$ID[0]\t$ID[1]\t$mutations{INH}\t$mutations{INH_Freq}\t$mutations{RMP}\t$mutations{RMP_Freq}\t$mutations{SM}\t$mutations{SM_Freq}\t$mutations{EMB}\t$mutations{EMB_Freq}\t$mutations{PZA}\t$mutations{PZA_Freq}\t$mutations{MFX}\t$mutations{MFX_Freq}\t$mutations{LFX}\t$mutations{LFX_Freq}\t$mutations{CFZ}\t$mutations{CFZ_Freq}\t$mutations{KAN}\t$mutations{KAN_Freq}\t$mutations{AMK}\t$mutations{AMK_Freq}\t$mutations{CPR}\t$mutations{CPR_Freq}\t$mutations{ETH}\t$mutations{ETH_Freq}\t$mutations{LZD}\t$mutations{LZD_Freq}\t$mutations{BDQ}\t$mutations{BDQ_Freq}\t$mutations{CS}\t$mutations{CS_Freq}\t$mutations{PAS}\t$mutations{PAS_Freq}\t$mutations{DLM}\t$mutations{DLM_Freq}\t$prediction\n";
    $R->stopR();
    close Fin;
    close Fout;
 }
}


sub tbcombinedresi{
	# get parameter and input from front-end.
   my $logprint                   =      shift;
   my $RESI_OUT                   =      shift;
   my $resi_list_date             =      shift;
   my @resisum_files              =      @_;
   my $resisum_file               =      "";
   my $id                         =      "";
   my @sample                     ;
   my $output_file                =  "Strain_Resistance.tab";
   my %check_up;
	
	#opendir(RESIDIR,"$RESI_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $RESI_OUT: MTBseq.pl line: ", __LINE__ ," \n";
	#@resisum_files       =  grep { $_ =~ /^\w.*\_${resi_list_date}\_resi\_summary.tab$/ && -f "$RESI_OUT/$_"   }  readdir(RESIDIR);

	
	if(-f "$RESI_OUT/$output_file") {
      open(IN,"$RESI_OUT/$output_file") || die print $logprint "\n<ERROR>\t",timer(),"\tCan't open $output_file: TBresi.pm line: ", __LINE__ , " \n";
      <IN>;
      while(<IN>) {
         my $line       =  $_;
         $line          =~ s/\015?\012?$//;
         my @fields     =  split(/\t/);
         $check_up{$fields[0]."_".$fields[1]}   =  1;
      }
	
	 close IN;
   }
  unless(-f "$RESI_OUT/$output_file") {
      print $logprint "<INFO>\t",timer(),"\t","Start writing $output_file...\n";
      open (Fout,">>$RESI_OUT/${output_file}") or die "\n<ERROR>\t",timer(),"Cannot open/write $output_file file\n";
      print Fout "SampleID\tLibID\tINH\tFreq_INH\tRMP\tFreq_RMP\tSM\tFreq_SM\tEMB\tFreq_EMB\tPZA\tFreq_PZA\tMFX\tFreq_MFX\tLFX\tFreq_LFX\tCFZ\tFreq_CFZ\tKAN\tFreq_KAN\tAMK\tFreq_AMK\tCPR\tFreq_CPR\tETH/PTH\tFreq_ETH/PTH\tLZD\tFreq_LZD\tBDQ\tFreq_BDQ\tCS\tFreq_CS\tPAS\tFreq_PAS\tDLM\tFreq_DLM\tPrediction\n";
      close Fout;
      print $logprint "<INFO>\t",timer(),"\t","Finished writing $output_file!\n";
   }
   open (Fout,">>$RESI_OUT/${output_file}") or die "\n<ERROR>\t",timer(),"Cannot open/write $output_file file: TBresi.pm line: ", __LINE__ , " \n";
   
   foreach my $resisum_file (sort { $a cmp $b } @resisum_files) {
	print $logprint "<INFO>\t",timer(),"\t","Start parsing $resisum_file...\n";
    $resisum_file       =~ /(\S+)\_resi_summary\.tab$/;
    my $id      =  $1;
    my @sample  =  split(/_/,$id);
	   
       # check if resistance summary already exists.
    if(exists $check_up{"$sample[0]"."_"."$sample[1]"}) {
         print $logprint "<INFO>\t",timer(),"\t","Skipping $resisum_file. Resistance classification already existing!\n";
         next;
		}
	  
    
    open (Fin,"<$RESI_OUT/$resisum_file") or die "\n<ERROR>\t",timer(),"Cannot open $resisum_file: TBresi.pm line: ", __LINE__ , " \n";
   
	my $line = <Fin>; #remove header
	my $mutations =<Fin>; #store mutation line in a variable
	print Fout "$mutations";
	close Fin;
	print $logprint "<INFO>\t",timer(),"\t","Finished writing resistance result for $id!\n";
    $mutations         =  {};
	$id                =  {};
	}
	
	undef(%check_up);
	close Fout;
}
1;
