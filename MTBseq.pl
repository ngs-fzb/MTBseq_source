#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use Cwd;
use Getopt::Long;
use TBbwa;
use TBstats;
use TBrefine;
use TBpile;
use TBlist;
use TBvariants;
use TBjoin;
use TBamend;
use TBstrains;
use TBgroups;
use TBtools;

# get current working directory and time.
my $W_dir         =     getcwd();
my $date_string   =     timer();
$date_string      =~    s/\[//;
$date_string      =~    s/\]//;
$date_string      =~    s/\-/_/g;
$date_string      =~    s/\:/_/g;
$date_string      =~    s/\s.*//;
$date_string      =~    s/_/-/g;

# define working directories.
my $BAM_OUT       =     "$W_dir/Bam";
my $GATK_OUT      =     "$W_dir/GATK_Bam";
my $MPILE_OUT     =     "$W_dir/Mpileup";
my $POS_OUT       =     "$W_dir/Position_Tables";
my $CALL_OUT      =     "$W_dir/Called";
my $STATS_OUT     =     "$W_dir/Statistics";
my $JOIN_OUT      =     "$W_dir/Joint";
my $AMEND_OUT     =     "$W_dir/Amend";
my $STRAIN_OUT    =     "$W_dir/Classification";
my $GROUPS_OUT    =     "$W_dir/Groups";

# define PATHs.
my $VAR_dir       =     "$RealBin/var/ref";
my $BWA_dir       =     "$RealBin/opt/bwa_0.7.17";
my $SAMTOOLS_dir  =     "$RealBin/opt/samtools_1.6";
my $PICARD_dir    =     "$RealBin/opt/picard_2.17.0";
my $IGV_dir       =     "$RealBin/opt/IGVTools_2.3.98";
my $GATK_dir      =     "$RealBin/opt/GenomeAnalysisTK_3.8";

# initialize command-line parameter.
my $step                =     "";
my $continue            =     "";
my $samples             =     "";
my $group_name          =     "";
my $resi_list_master    =     "";
my $int_regions         =     "";
my $categories          =     "";
my $basecalib           =     "";
my $ref                 =     "";
my $mibqual             =     "";
my $all_vars            =     "";
my $snp_vars            =     "";
my $lowfreq_vars        =     "";
my $micovf              =     "";
my $micovr              =     "";
my $miphred20           =     "";
my $mifreq              =     "";
my $unambigous          =     "";
my $window              =     "";
my $distance            =     "";
my $quiet               =     "";
my $threads             =     "";
my $naming_scheme       =     "";
my $help                =     "";

# get command-line parameter.
GetOptions('step:s'        =>    \$step,
           'continue'      =>    \$continue,
           'samples:s'     =>    \$samples,
           'project:s'     =>    \$group_name,
           'resilist:s'    =>    \$resi_list_master,
           'intregions:s'  =>    \$int_regions,
           'categories:s'  =>    \$categories,
           'basecalib:s'   =>    \$basecalib,
           'ref:s'         =>    \$ref,
           'minbqual:i'    =>    \$mibqual,
           'all_vars'      =>    \$all_vars,
           'snp_vars'      =>    \$snp_vars,
           'lowfreq_vars'  =>    \$lowfreq_vars,
           'mincovf:i'     =>    \$micovf,
           'mincovr:i'     =>    \$micovr,
           'minphred20:i'  =>    \$miphred20,
           'minfreq:i'     =>    \$mifreq,
           'unambig:i'     =>    \$unambigous,
           'window:i'      =>    \$window,
           'distance:i'    =>    \$distance,
           'quiet'         =>    \$quiet,
           'threads:i'     =>    \$threads,
           'help'          =>    \$help
          );

# print help message if specified or error if step is not defined or wrong.
if($help eq '1') { help();     exit 1; }
if($step eq '' ) { nostep();   exit 1; }
unless(($step eq 'TBfull'     )  ||
       ($step eq 'TBbwa'      )  ||
       ($step eq 'TBstats'    )  ||
       ($step eq 'TBrefine'   )  ||
       ($step eq 'TBpile'     )  ||
       ($step eq 'TBlist'     )  ||
       ($step eq 'TBvariants' )  ||
       ($step eq 'TBjoin'     )  ||
       ($step eq 'TBamend'    )  ||
       ($step eq 'TBstrains'  )  ||
       ($step eq 'TBgroups'   )
      ) { badstep($step); exit 1; }

# test command-line options and set default vaules if necessary.
if($continue         eq    ''    ) { $continue           =     0;                                  }
if($samples          eq    ''    ) { $samples            =     "NONE";                             }
if($group_name       eq    ''    ) { $group_name         =     "NONE";                             }
if($resi_list_master eq    ''    ) { $resi_list_master   =     "NONE";                             }
if($int_regions      eq    ''    ) { $int_regions        =     "NONE";                             }
if($categories       eq    ''    ) { $categories         =     "NONE";                             }
if($basecalib        eq    ''    ) { $basecalib          =     "NONE";                             }
if($ref              eq    ''    ) { $ref                =     "M._tuberculosis_H37Rv_2015-11-13"; }
if($mibqual          eq    ''    ) { $mibqual            =     13;                                 }
if($all_vars         eq    ''    ) { $all_vars           =     0;                                  }
if($snp_vars         eq    ''    ) { $snp_vars           =     0;                                  }
if($lowfreq_vars     eq    ''    ) { $lowfreq_vars       =     0;                                  }
if($micovf           eq    ''    ) { $micovf             =     4;                                  }
if($micovr           eq    ''    ) { $micovr             =     4;                                  }
if($miphred20        eq    ''    ) { $miphred20          =     4;                                  }
if($mifreq           eq    ''    ) { $mifreq             =     75;                                 }
if($unambigous       eq    ''    ) { $unambigous         =     95;                                 }
if($window           eq    ''    ) { $window             =     12;                                 }
if($distance         eq    ''    ) { $distance           =     12;                                 }
if($quiet            eq    ''    ) { $quiet              =     0;                                  }
if($threads          eq    ''    ) { $threads            =     1;                                  }
if($threads          >     8     ) { $threads            =     8;                                  }

# set name of $ref fasta file and gene annotation.
my $refg    =   $ref;
$ref        .=  ".fasta";
$refg       .=  "_genes.txt";

# if $ref is MTB than turn on resistance check.
if($ref eq 'M._tuberculosis_H37Rv_2015-11-13.fasta') {
   $resi_list_master    =   "$RealBin/var/res/MTB_Resistance_Mediating.txt"             if($resi_list_master eq 'NONE');   # MTBC Resistance mediating variants and MTBC phylogentic SNPs.
   $int_regions         =   "$RealBin/var/res/MTB_Extended_Resistance_Mediating.txt"    if($int_regions eq 'NONE');        # MTBC extended intergenic resistance mediating positions.
   $categories          =   "$RealBin/var/cat/MTB_Gene_Categories.txt"                  if($categories eq 'NONE');         # MTBC essential and non-essential genes.
   $basecalib           =   "$RealBin/var/res/Base_Calibration_List.vcf"                if($basecalib eq 'NONE');          # MTBC base calibration list. Known SNP positions for base call recalibration.
}

# create log file and log on screen if --quiet is unset.
my $logprint;
if($quiet == 0) {
   open($logprint, "|-", "tee -a $W_dir/MTBseq_$date_string\_$ENV{USER}.log");
}

# and if --quiet is set.
else {
   open($logprint, ">>", "$W_dir/MTBseq_$date_string\_$ENV{USER}.log");
}
select($logprint);
$| = 1;

# inform the user what you will do.
print $logprint "\n<INFO>\t",timer(),"\tYou are $ENV{USER}.\n";
print $logprint "<INFO>\t",timer(),"\tYour current working directory is: $W_dir\n";
print $logprint "<INFO>\t",timer(),"\tYou requested $threads thread(s) for the pipeline.\n";
print $logprint "\n<INFO>\t",timer(),"\tYour parameter setting is:\n";
print $logprint "<INFO>\t",timer(),"\t--step\t\t$step\n";
print $logprint "<INFO>\t",timer(),"\t--continue\t$continue\n";
print $logprint "<INFO>\t",timer(),"\t--samples\t$samples\n";
print $logprint "<INFO>\t",timer(),"\t--project\t$group_name\n";
print $logprint "<INFO>\t",timer(),"\t--resilist\t$resi_list_master\n";
print $logprint "<INFO>\t",timer(),"\t--intregions\t$int_regions\n";
print $logprint "<INFO>\t",timer(),"\t--categories\t$categories\n";
print $logprint "<INFO>\t",timer(),"\t--basecalib\t$basecalib\n";
print $logprint "<INFO>\t",timer(),"\t--ref\t\t$ref\n";
print $logprint "<INFO>\t",timer(),"\t--minbqual\t$mibqual\n";
print $logprint "<INFO>\t",timer(),"\t--all_vars\t$all_vars\n";
print $logprint "<INFO>\t",timer(),"\t--snp_vars\t$snp_vars\n";
print $logprint "<INFO>\t",timer(),"\t--lowfreq_vars\t$lowfreq_vars\n";
print $logprint "<INFO>\t",timer(),"\t--mincovf\t$micovf\n";
print $logprint "<INFO>\t",timer(),"\t--mincovr\t$micovr\n";
print $logprint "<INFO>\t",timer(),"\t--minphred20\t$miphred20\n";
print $logprint "<INFO>\t",timer(),"\t--minfreq\t$mifreq\n";
print $logprint "<INFO>\t",timer(),"\t--unambig\t$unambigous\n";
print $logprint "<INFO>\t",timer(),"\t--window\t$window\n";
print $logprint "<INFO>\t",timer(),"\t--distance\t$distance\n";
print $logprint "<INFO>\t",timer(),"\t--quiet\t\t$quiet\n";

print $logprint "\n<INFO>\t",timer(),"\tThe following programs will be used, if necessary:\n"; 
print $logprint "<INFO>\t",timer(),"\t$BWA_dir\n";
print $logprint "<INFO>\t",timer(),"\t$SAMTOOLS_dir\n";
print $logprint "<INFO>\t",timer(),"\t$PICARD_dir\n";
print $logprint "<INFO>\t",timer(),"\t$IGV_dir\n";
print $logprint "<INFO>\t",timer(),"\t$GATK_dir\n";

print $logprint "\n<INFO>\t",timer(),"\tThe following directories will be used, if necessary:\n";
print $logprint "<INFO>\t",timer(),"\t$BAM_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$GATK_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$MPILE_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$POS_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$CALL_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$JOIN_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$AMEND_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$STRAIN_OUT\n";
print $logprint "<INFO>\t",timer(),"\t$GROUPS_OUT\n";

system("mkdir -p $BAM_OUT");
system("mkdir -p $GATK_OUT");
system("mkdir -p $MPILE_OUT");
system("mkdir -p $POS_OUT");
system("mkdir -p $CALL_OUT");
system("mkdir -p $STATS_OUT");
system("mkdir -p $JOIN_OUT");
system("mkdir -p $AMEND_OUT");
system("mkdir -p $STRAIN_OUT");
system("mkdir -p $GROUPS_OUT");

# Initialize check up and content arrays.
my %check_up;
my @fastq_files;
my @fastq_files_new;
my @bam_files;
my @bam_files_new;
my @gatk_files;
my @gatk_files_new;
my @mpile_files;
my @mpile_files_new;
my @pos_files;
my @pos_files_new;
my @var_files;
my @var_files_new;
my @join_files;
my @join_files_new;
my @amend_files;
my @amend_files_new;
my @group_files;
my @samples;
my $sample_number;
my $output_mode = $all_vars . $snp_vars . $lowfreq_vars;
if(-f $samples) {
   open(IN,"$samples") || die print $logprint "<INFO>\t",timer(),"\tCan't find $samples file! MTBseq.pl line: ", __LINE__ ," \n";
   while(my $line = <IN>) {
      chomp($line);
      $line =~ s/\015?\012?$//;
      next unless ($line);
      $sample_number++;
   }
   close(IN);
}

# jump to certain pipeline steps.
if($step eq 'TBfull'       ) { goto TBfull      }
if($step eq 'TBbwa'        ) { goto TBbwa       }
if($step eq 'TBrefine'     ) { goto TBrefine    }
if($step eq 'TBpile'       ) { goto TBpile      }
if($step eq 'TBlist'       ) { goto TBlist      }
if($step eq 'TBvariants'   ) { goto TBvariants  }
if($step eq 'TBstats'      ) { goto TBstats     }
if($step eq 'TBjoin'       ) { goto TBjoin      }
if($step eq 'TBamend'      ) { goto TBamend     }
if($step eq 'TBstrains'    ) { goto TBstrains   }
if($step eq 'TBgroups'     ) { goto TBgroups    }

# Pipeline execution.
TBfull:
print $logprint "\n<INFO>\t",timer(),"\t### [TBfull] selected ###\n";



TBbwa:
if($step eq 'TBbwa') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBbwa] selected ###\n";
}
opendir(WORKDIR,"$W_dir")       || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $W_dir: MTBseq.pl line: ", __LINE__ ," \n";
opendir(BAMDIR,"$BAM_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $BAM_OUT: MTBseq.pl line: ", __LINE__ , "\n";
@fastq_files      =  grep { $_ =~ /^\w.*R\d+\.fastq\.gz/ && -f "$W_dir/$_"    }  readdir(WORKDIR);
@bam_files        =  grep { $_ =~ /^\w.*\.bam$/ && -f "$BAM_OUT/$_"           }  readdir(BAMDIR);
closedir(WORKDIR);
closedir(BAMDIR);
if(scalar(@fastq_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo read files to map! Check content of $W_dir!\n";
   exit 1;
}
%check_up         =  map { (my $id = $_) =~ s/\.bam$//; $id => $id; } @bam_files;
for(my $i = 0; $i < scalar(@fastq_files); $i++) {
   my $tmp        =  $fastq_files[$i];
   $tmp           =~ s/_R\d\.fastq.gz//;
   if(!exists $check_up{$tmp}) {
      push(@fastq_files_new,$fastq_files[$i]);
   }
}
if(scalar(@fastq_files_new) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tAll read files have been already mapped from $W_dir!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tMapping samples:\n";
foreach my $fastq (sort { $a cmp $b } @fastq_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$fastq\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart BWA mapping...\n";
tbbwa($logprint,$W_dir,$VAR_dir,$BWA_dir,$SAMTOOLS_dir,$BAM_OUT,$ref,$threads,$naming_scheme,@fastq_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished BWA mapping!\n";
@fastq_files      =  ();
@bam_files        =  ();
@fastq_files_new  =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBrefine:
if($step eq 'TBrefine') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBrefine] selected ###\n";
}
opendir(BAMDIR,"$BAM_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $BAM_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(GATKDIR,"$GATK_OUT")    || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $GATK_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@bam_files        =  grep { $_ =~ /^\w.*\.bam$/ && -f "$BAM_OUT/$_"        }  readdir(BAMDIR);
@gatk_files       =  grep { $_ =~ /^\w.*\.gatk\.bam$/ && -f "$GATK_OUT/$_" }  readdir(GATKDIR);
closedir(BAMDIR);
closedir(GATKDIR);
if(scalar(@bam_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo mapping files to refine! Check content of $BAM_OUT!\n";
   exit 1;
}
%check_up         =  map { (my $id = $_) =~ s/\.gatk\.bam$//; $id => $id; } @gatk_files;
for(my $i = 0; $i < scalar(@bam_files); $i++) {
   my $tmp        =  $bam_files[$i];
   $tmp           =~ s/\.bam//;
   if(!exists $check_up{$tmp}) {
      push(@gatk_files_new,$bam_files[$i]);
   }
}
if(scalar(@gatk_files_new) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tAll mapping files have been already refined from $BAM_OUT!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tRefining mappings:\n";
foreach my $bam (sort { $a cmp $b } @gatk_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$bam\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart GATK refinement...\n";
tbrefine($logprint,$W_dir,$VAR_dir,$PICARD_dir,$IGV_dir,$GATK_dir,$BAM_OUT,$GATK_OUT,$ref,$basecalib,$threads,@gatk_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished GATK logic!\n";
@bam_files        =  ();
@gatk_files       =  ();
@gatk_files_new   =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBpile:
if($step eq 'TBpile') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBpile] selected ###\n";
}
opendir(GATKDIR,"$GATK_OUT")    || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $GATK_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(MPILEDIR,"$MPILE_OUT")  || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $MPILE_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@gatk_files       =  grep { $_ =~ /^\w.*\.gatk\.bam$/ && -f "$GATK_OUT/$_"       }  readdir(GATKDIR);
@mpile_files      =  grep { $_ =~ /^\w.*\.gatk\.mpileup$/ && -f "$MPILE_OUT/$_"  }  readdir(MPILEDIR);
closedir(GATKDIR);
closedir(MPILEDIR);
if(scalar(@gatk_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo files to create pileups! Check content of $GATK_OUT!\n";
   exit 1;
}
%check_up         =  map { (my $id = $_) =~ s/\.mpileup$//; $id => $id; } @mpile_files;
for(my $i = 0; $i < scalar(@gatk_files); $i++) {
   my $tmp        =  $gatk_files[$i];
   $tmp           =~ s/\.bam//;
   if(!exists $check_up{$tmp}) {
      push(@gatk_files_new,$gatk_files[$i]);
   }
}
if(scalar(@gatk_files_new) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tAll mpileup files have been already created from $GATK_OUT!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tCreating mpileups:\n";
foreach my $bam (sort { $a cmp $b } @gatk_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$bam\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart creating mpileup files...\n";
tbpile($logprint,$VAR_dir,$SAMTOOLS_dir,$GATK_dir,$GATK_OUT,$MPILE_OUT,$ref,$threads,@gatk_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished creating mpileup files!\n";
@gatk_files       =  ();
@mpile_files      =  ();
@gatk_files_new   =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBlist:
if($step eq 'TBlist') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBlist] selected ###\n";
}
opendir(MPILEDIR,"$MPILE_OUT")  || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $MPILE_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(POSDIR,"$POS_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $POS_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@mpile_files      =  grep { $_ =~ /^\w.*\.gatk\.mpileup$/ && -f "$MPILE_OUT/$_"           }   readdir(MPILEDIR);
@pos_files        =  grep { $_ =~ /^\w.*\.gatk_position_table\.tab$/ && -f "$POS_OUT/$_"  }   readdir(POSDIR);
closedir(MPILEDIR);
closedir(POSDIR);
if(scalar(@mpile_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo files to create position lists! Check content of $MPILE_OUT!\n";
   exit 1;
}
%check_up         =  map { (my $id = $_) =~ s/\.gatk_position_table\.tab$//; $id => $id; } @pos_files;
for(my $i = 0; $i < scalar(@mpile_files); $i++) {
   my $tmp        =  $mpile_files[$i];
   $tmp           =~ s/\.gatk\.mpileup//;
   if(!exists $check_up{$tmp}) {
      push(@mpile_files_new,$mpile_files[$i]);
   }
}
if(scalar(@mpile_files_new) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tAll position lists have been already created from $MPILE_OUT!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tCreating position lists:\n";
foreach my $mpile (sort { $a cmp $b } @mpile_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$mpile\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart creating position lists...\n";
tblist($logprint,$VAR_dir,$MPILE_OUT,$POS_OUT,$ref,$mibqual,$threads,@mpile_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished creating position lists!\n";
@mpile_files      =  ();
@pos_files        =  ();
@mpile_files_new  =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBvariants:
if($step eq 'TBvariants') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBvariants] selected ###\n";
}
opendir(POSDIR,"$POS_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $POS_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(CALLDIR,"$CALL_OUT")    || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $CALL_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@pos_files        =  grep { $_ =~ /^\w.*\.gatk_position_table\.tab$/ && -f "$POS_OUT/$_"                                                                          }  readdir(POSDIR);
@var_files        =  grep { $_ =~ /^\w.*\.gatk_position_variants_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode$output_mode\.tab$/ && -f "$CALL_OUT/$_"   }  readdir(CALLDIR);
closedir(POSDIR);
closedir(CALLDIR);
if(scalar(@pos_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo files to call variants! Check content of $POS_OUT!\n";
   exit 1;
}
%check_up         =  map { (my $id = $_) =~ s/\.gatk_position_variants_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_outmode$output_mode\.tab$//; $id => $id; } @var_files;
for(my $i = 0; $i < scalar(@pos_files); $i++) {
   my $tmp        =  $pos_files[$i];
   $tmp           =~ s/\.gatk_position_table\.tab//;
   if(!exists $check_up{$tmp}) {
      push(@pos_files_new,$pos_files[$i]);
   }
}
if(scalar(@pos_files_new) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tAll variant files have been already created from $POS_OUT!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tCalling variants:\n";
foreach my $pos_file (sort { $a cmp $b } @pos_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$pos_file\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart variant calling...\n";
tbvariants($logprint,$VAR_dir,$POS_OUT,$CALL_OUT,$ref,$refg,$micovf,$micovr,$miphred20,$mifreq,$resi_list_master,$int_regions,$all_vars,$snp_vars,$lowfreq_vars,@pos_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished variant calling!\n";
@pos_files     =  ();
@var_files     =  ();
@pos_files_new =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBstats:
if($step eq 'TBstats') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBstats] selected ###\n";
}
opendir(BAMDIR,"$BAM_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $BAM_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(POSDIR,"$POS_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $POS_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@bam_files     =  grep { $_ =~ /^\w.*\.bam$/ && -f "$BAM_OUT/$_"                          }  readdir(BAMDIR);
@pos_files     =  grep { $_ =~ /^\w.*\.gatk_position_table\.tab$/ && -f "$POS_OUT/$_"     }  readdir(POSDIR);
closedir(BAMDIR);
closedir(POSDIR);
if(scalar(@bam_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo files to calculate statistics! Check content of $BAM_OUT!\n";
   exit 1;
}
%check_up      =  map { (my $id = $_) =~ s/\.gatk_position_table\.tab$//; $id => $id; } @pos_files;
for(my $i = 0; $i < scalar(@bam_files); $i++) {
   my $tmp     =  $bam_files[$i];
   $tmp        =~ s/\.bam//;
   if(exists $check_up{$tmp}) {
      push(@bam_files_new,$bam_files[$i]);
   }
   else {
      print $logprint "\n<INFO>\t",timer(),"\tThe following mapping has no position list: $bam_files[$i]. This will skip the statistics part for this sample!\n";
   }
}
print $logprint "\n<INFO>\t",timer(),"\tCalculating statistics for mappings:\n";
foreach my $bam (sort { $a cmp $b } @bam_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$bam\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart statistics calculation...\n";
tbstats($logprint,$W_dir,$VAR_dir,$SAMTOOLS_dir,$BAM_OUT,$POS_OUT,$STATS_OUT,$ref,$refg,$micovf,$micovr,$miphred20,$mifreq,$all_vars,$snp_vars,$lowfreq_vars,$date_string,@bam_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished statistics calculation!\n";
@bam_files     =  ();
@pos_files     =  ();
@bam_files_new =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBstrains:
if($step eq 'TBstrains') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBstrains] selected ###\n";
}
opendir(POSDIR,"$POS_OUT")      || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $POS_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@pos_files     =  grep { $_ =~ /^\w.*\.gatk_position_table\.tab$/ && -f "$POS_OUT/$_" } readdir(POSDIR);
closedir(POSDIR);
if(scalar(@pos_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo files to call strains! Check content of $POS_OUT!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tCalling strains:\n";
foreach my $pos_file (sort { $a cmp $b } @pos_files) {
   print $logprint "<INFO>\t",timer(),"\t$pos_file\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart strain call...\n";
tbstrains($logprint,$VAR_dir,$POS_OUT,$STRAIN_OUT,$ref,$refg,$micovf,$micovr,$mifreq,$miphred20,$date_string,$all_vars,$snp_vars,$lowfreq_vars,@pos_files);
print $logprint "<INFO>\t",timer(),"\tFinished strain call!\n";
@pos_files     =  ();
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBjoin:
if($step eq 'TBjoin') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBjoin] selected ###\n";
}
if($samples eq 'NONE') {
   print $logprint "\n<INFO>\t",timer(),"\tSkipping Joint analysis. No samples defined for joint analysis.\n";
   print $logprint "\n<INFO>\t",timer(),"\t### MTBseq finished!!! ###\n";
   exit 1;
}
open(IN,"$samples") || die print $logprint "<INFO>\t",timer(),"\tCan't find $samples file!\n";
while(<IN>) {
   my $line    =  $_;
   $line       =~ s/\015?\012?$//;
   next unless ($line);
   my @fields  =  split(/\t/,$line);
   $check_up{$fields[0]."_".$fields[1]}   =  $fields[0]."_".$fields[1];
}
close(IN);
opendir(CALLDIR,"$CALL_OUT")    || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $CALL_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(JOINDIR,"$JOIN_OUT")    || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $JOIN_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@var_files     =  grep { $_ =~ /^\w.*\.gatk_position_variants_cf$micovf\_cr$micovf\_fr$mifreq\_ph$miphred20\_outmode$output_mode\.tab$/ && -f "$CALL_OUT/$_"   }  readdir(CALLDIR);
@join_files    =  grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_samples$sample_number\.tab$/ && -f "$JOIN_OUT/$_"            }  readdir(JOINDIR);
closedir(CALLDIR);
closedir(JOINDIR);
if(scalar(@var_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo files to create joint variant tables! Check content of $POS_OUT and $CALL_OUT!\n";
	exit 1;
}
if(scalar(@join_files) != 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tIt seems that you already created a joint analysis for project $group_name! Check content of $JOIN_OUT!\n";
   exit 1;
}
for(my $i = 0; $i < scalar(@var_files); $i++) {
   my @tmp     =  split(/_|\./,$var_files[$i]);
   if(exists $check_up{$tmp[0]."_".$tmp[1]}) {
      print $logprint "\n<INFO>\t",timer(),"\tFound $var_files[$i]! Using file for joint analysis...";
      push(@var_files_new,$var_files[$i]);
      delete($check_up{$tmp[0]."_".$tmp[1]});
   }
}
foreach my $samplelib(keys %check_up) {
   print $logprint "\n<ERROR>\t",timer(),"\t Could not find files for $samplelib! Did you run TBlist and TBvariants for $samplelib?\n";
   exit 1;
}
print $logprint "\n\n<INFO>\t",timer(),"\tCreating a joint variant table from:\n";
foreach my $var_file (sort { $a cmp $b } @var_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$var_file\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart creating joint variant table...\n";
tbjoin($logprint,$VAR_dir,$POS_OUT,$CALL_OUT,$JOIN_OUT,$group_name,$ref,$refg,$micovf,$micovr,$miphred20,$mifreq,$all_vars,$snp_vars,$lowfreq_vars,@var_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished creating joint variant table!\n";
@var_files     =  ();
@join_files    =  ();
@var_files_new =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBamend:
if($step eq 'TBamend') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBamend] selected ###\n";
}
opendir(JOINDIR,"$JOIN_OUT")    || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $JOIN_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(AMENDDIR,"$AMEND_OUT")  || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $AMEND_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@join_files    =  grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_samples$sample_number\.tab$/ && -f "$JOIN_OUT/$_"                                       }  readdir(JOINDIR);
@amend_files   =  grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_samples$sample_number\_amended_u$unambigous\_phylo_w$window.tab$/ && -f "$AMEND_OUT/$_" }  readdir(AMENDDIR);
closedir(JOINDIR);
closedir(AMENDDIR);
if(scalar(@join_files) == 0) {
   my $joint_file_name_tmp = $group_name . "_joint_cf" . $micovf . "_cr" . $micovr . "_fr" . $mifreq . "_ph" . $miphred20 . "_samples" . $sample_number . ".tab";
   print $logprint "\n<ERROR>\t",timer(),"\tNo joint variant file $joint_file_name_tmp to amend! Check content of $JOIN_OUT!\n";
   exit 1;
}
%check_up      =  map { (my $id = $_) =~ s/\_amended_u$unambigous\_phylo_w$window.tab//; $id => $id; } @amend_files;
for(my $i = 0; $i < scalar(@join_files); $i++) {
   my $tmp     =  $join_files[$i];
   $tmp        =~ s/\.tab//;
   if(!exists $check_up{$tmp}) {
      push(@join_files_new,$join_files[$i]);
   }
}
if(scalar(@join_files_new) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tA Joint variant file with project name \"$group_name\" is already existing! Check content of $JOIN_OUT!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tAmending joint variant table:\n";
foreach my $join_file (sort { $a cmp $b } @join_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$join_file\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart amending joint variant table...\n";
tbamend($logprint,$JOIN_OUT,$AMEND_OUT,$micovf,$micovr,$miphred20,$mifreq,$unambigous,$window,$categories,$resi_list_master,$int_regions,@join_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished amending joint variant table!\n";
@join_files       =  ();
@amend_files      =  ();
@join_files_new   =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



TBgroups:
if($step eq 'TBgroups') {
   print $logprint "\n<INFO>\t",timer(),"\t### [TBgroups] selected ###\n";
}
opendir(AMENDDIR,"$AMEND_OUT")  || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $AMEND_OUT: MTBseq.pl line: ", __LINE__ ," \n";
opendir(GROUPDIR,"$GROUPS_OUT") || die print $logprint "<ERROR>\t",timer(),"\tCan\'t open directory $GROUPS_OUT: MTBseq.pl line: ", __LINE__ ," \n";
@amend_files      =  grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_samples$sample_number\_amended_u$unambigous\_phylo_w$window\.tab$/ && -f "$AMEND_OUT/$_"      }readdir(AMENDDIR);
@group_files      =  grep { $_ =~ /$group_name\_joint_cf$micovf\_cr$micovr\_fr$mifreq\_ph$miphred20\_samples$sample_number\_amended_u$unambigous\_phylo_w$window\.matrix$/ && -f "$GROUPS_OUT/$_"  }  readdir(GROUPDIR);
closedir(AMENDDIR);
closedir(GROUPDIR);
if(scalar(@amend_files) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tNo file to call groups! Check content of $AMEND_OUT and the project name $group_name!\n";
   exit 1;
}
%check_up         =  map { (my $id = $_) =~ s/\.matrix$//; $id => $id; } @group_files;
for(my $i = 0; $i < scalar(@amend_files); $i++) {
   my $tmp        =  $amend_files[$i];
   $tmp           =~ s/\.tab//;
   if(!exists $check_up{$tmp}) {
      push(@amend_files_new,$amend_files[$i]);
   }
}
if(scalar(@amend_files_new) == 0) {
   print $logprint "\n<ERROR>\t",timer(),"\tA group call with project name \"$group_name\" is already existing! Check content of $GROUPS_OUT!\n";
   exit 1;
}
print $logprint "\n<INFO>\t",timer(),"\tCalling groups from:\n";
foreach my $amend_file (sort { $a cmp $b } @amend_files_new) {
   print $logprint "<INFO>\t",timer(),"\t$amend_file\n";
}
print $logprint "\n<INFO>\t",timer(),"\tStart group calling...\n";
tbgroups($logprint,$AMEND_OUT,$GROUPS_OUT,$distance,@amend_files_new);
print $logprint "<INFO>\t",timer(),"\tFinished group calling!\n";
@amend_files      =  ();
@group_files      =  ();
@amend_files_new  =  ();
undef(%check_up);
if($continue == 0 && $step ne 'TBfull') { exit 1; }



print $logprint "\n<INFO>\t",timer(),"\t### MTBseq finished!!! ###\n";
close($logprint);
