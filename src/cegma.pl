#!/usr/bin/perl

##########################################################################
#                                                                        #
#                               CEGMA                                    #
#                                                                        #
##########################################################################
#                                                                        #
#                Core Eukaryotic Gene Mapping Approach.                  #
#                                                                        #
#              Copyright (C) 2006-2014       Genis Parra                 #
#                                            Keith Bradnam               #
#                                            Ian Korf                    #
#                                                                        #
# This program is free software; you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation; either version 2 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# This program is distributed in the hope that it will be useful, but    #
# WITHOUT ANY WARRANTY; without even the implied warranty of             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      #
# General Public License for more details.                               #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with this program; if not, write to the Free Software            #
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#                                                                        #
##########################################################################

use strict; use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Cegma;

##########################################################################
##                                                                      ##
##                     VARIABLES  AND DEFINITIONS                       ##
##                                                                      ##
##########################################################################

# MAIN VARIABLES 
my $PROGRAM = "cegma";
my $VERSION = "2.5";
 
my $CEGMA     = defined($ENV{'CEGMA'})    ?  $ENV{'CEGMA'} : '.';
my $BIN       = defined($ENV{'CEGMA'})    ? "$CEGMA/bin" : '.';
my $CEGMAdata = defined($ENV{'CEGMA'})    ? "$CEGMA/data" : "../data";
my $TMP       = defined($ENV{'CEGMATMP'}) ?  $ENV{'CEGMATMP'} : '/tmp';
my $TMPROOT   = "cegma_$$";
my $CEGMATMP  = "$TMP/$TMPROOT";


##########################################################################
##                                                                      ##
##                          READING OPTIONS                             ##
##                                                                      ##
##########################################################################

# print help and exit if no command-line options are given
&print_help unless @ARGV; 

# getOptions Variables
my ( $genome, $protein, $output, $geneid_opt, $geneid_param, $blast_opt, $score_cutoff, 
	 $genewise_file, $tblastn_gff, $blastdb, $annotation, $ofn, $ps_output, $max_intron, 
	 $hmm_profiles, $hmm_prefix, $prot_num, $cutoff_file, $verbose_str,
	 $savefiles_flg, $verbose_flg, $ext_flg, $mam_flg, $vrt_flg, $quiet_flg, $temp_flg, 
	 $help_flg, $interlen, $boundaries, $threads, 
     ) = ( undef, undef, undef, undef, undef, undef, undef,
           undef, undef, undef, undef, undef, undef, undef,
           undef, undef, undef, undef, undef,
           0, 0, 0, 0, 0, 0, 0, 0, 5000, 2000, 1 );


# reading options
&Which_Options();


##########################################################################
##                                                                      ##
##                          INITIAL CHECKS                              ##
##                                                                      ##
##########################################################################

# Cheking if the external programs are in the path.
foreach my $tool (qw(makeblastdb tblastn genewise geneid hmmsearch)){
	&go_to_die() if not Cegma::check_tool_exists($tool);
}



##########################################################################
##                                                                      ##
##                          COMPUTING DATA                              ##
##                                                                      ##
##########################################################################


##########################################################################
#       Running genome_map.pl (unless genome chunks file is present)
##########################################################################
&print_header("MAPPING PROTEINS TO GENOME (TBLASTN)");
&run_genome_map($genome, $protein, $tblastn_gff, $blastdb, $annotation, $interlen, $boundaries); 	



##########################################################################
# Running local_map.pl, can skip if geneid file exists, otherwise can also use existing genewise file if present
##########################################################################
&print_header("MAKING INITIAL GENE PREDICTIONS FOR CORE GENES (GENEWISE + GENEID)");
if (-s "local.geneid.fa"){
	my $geneid_count = Cegma::count_seqs_in_fasta_file("local.geneid.fa");
	warn "NOTE: local.geneid.fa file already exists and contains $geneid_count predictions. Will use this file instead and skip running local_map.\n" unless ($quiet_flg);
} elsif (defined $genewise_file){
	warn "NOTE: Will use specified $genewise_file when running local_map.\n" unless ($quiet_flg);
	&run_local_map("local", "-g $genewise_file -f");	
} else{
	&run_local_map("local", "-f");
}



##########################################################################
# Filtering the set of predicted proteins 
##########################################################################
&print_header("FILTERING INITIAL PROTEINS PRODUCED BY GENEID (HMMER)");
if (-s "local.geneid.selected.fa"){
	my $count = Cegma::count_seqs_in_fasta_file("local.geneid.selected.fa");
	warn "NOTE: local.geneid.selected.fa file already exists and contains $count predictions. Will use this file instead and skip running hmm_select.\n" unless ($quiet_flg);	
} else{
	my $options = "-i $hmm_prefix -o local -t $threads ";
	$options .= "-v " if ($verbose_flg);
    my $command = "hmm_select $options $hmm_profiles local.geneid.fa $cutoff_file 2>$output.cegma.errors";
    &print_command("$command");
	system("$BIN/$command") && die "Can't run hmm_select\n";
	my $count = Cegma::count_seqs_in_fasta_file("local.geneid.selected.fa");
	warn "NOTE: Found $count geneid predictions with scores above threshold value\n" unless ($quiet_flg);	

	# may as well stop the script if there are zero geneid predictions
	prematurely_end_program() if($count == 0);
}

prematurely_end_program() if (-z "local.geneid.selected.gff");



##########################################################################
# Calculating geneid parameters from selected gene predictions
##########################################################################
&print_header("CALCULATING GENEID PARAMETERS FROM SELECTED GENEID PREDICTIONS");
if (-s "geneid_params/self.param"){
	warn "NOTE: geneid_params/self.param file already exists. Will use this file instead and skip running geneid-train.\n" unless ($quiet_flg);	
} else{
	my $options = "";
	$options = "-v" if ($verbose_flg);
	my $command = "geneid-train $options local.geneid.selected.gff local.geneid.selected.dna geneid_params 2>$output.cegma.errors";
	&print_command("$command");
	system ("$BIN/$command")  && &go_to_die("geneid-train did not work properly");

	if ($max_intron) {
		open(INTRON, ">", "geneid_params/intron.max") or die "Can't write to geneid_params/intron.max\n";
	    print INTRON "First+:Internal+             Internal+:Terminal+          30:$max_intron block\n";
	    close(INTRON);
	}

	$command = "make_paramfile $CEGMA/data/self.param.template \\
	    geneid_params/coding.initial.5.logs geneid_params/coding.transition.5.logs \\
	    geneid_params/start.logs geneid_params/acc.logs geneid_params/don.logs \\
	    geneid_params/intron.max >  geneid_params/self.param";

	&print_command("$command");
	system ("$BIN/$command") && &go_to_die("make_paramfile did not work properly.");	
}



##########################################################################
# Remapping more accurately based on updated geneid paramaters???
##########################################################################
&print_header("ACCURATE LOCAL MAPPING");

if (-s "local_self.geneid.fa"){
	my $geneid_count = Cegma::count_seqs_in_fasta_file("local_self.geneid.fa");
	warn "NOTE: local_self.geneid.fa file already exists and contains $geneid_count predictions. Will use this file instead and skip running local_map.\n" unless ($quiet_flg);	
} else{
	&run_local_map("local_self","-g local.genewise.gff -d geneid_params/self.param");
}



##########################################################################
# Filtering the set of predicted proteins 
##########################################################################
&print_header("FINAL FILTERING");

if (-s "local_self.geneid.selected.fa"){
	my $count = Cegma::count_seqs_in_fasta_file("local_self.geneid.selected.fa");
	warn "NOTE: local_self.geneid.selected.fa file already exists and contains $count predictions. Will use this file instead and skip running hmm_select.\n" unless ($quiet_flg);	
} else{
	my $options = "-i $hmm_prefix -o local_self -t $threads ";
	$options .= "-v " if ($verbose_flg);
    my $command = "hmm_select $options $hmm_profiles local_self.geneid.fa $cutoff_file 2>$output.cegma.errors";

    &print_command("$command");
	system("$BIN/$command") && die "Can't run $command\n";
	my $count = Cegma::count_seqs_in_fasta_file("local_self.geneid.selected.fa");
	prematurely_end_program() if ($count == 0);
	warn "NOTE: Found $count geneid predictions with scores above threshold value\n" unless ($quiet_flg);	
}


# rename files final to output names
system("sed 's/geneid_v.../cegma/' < local_self.geneid.selected.gff > $output.cegma.local.gff") && die "Can't run sed\n";
rename("local_self.geneid.selected.id",  "$output.cegma.id")  or die "Can't rename file\n";
rename("local_self.geneid.selected.fa",  "$output.cegma.fa")  or die "Can't rename file\n";
rename("local_self.geneid.selected.dna", "$output.cegma.dna") or die "Can't rename file\n";



##########################################################################
# Convert local coordinates into genome wide coordinates
##########################################################################
&print_header("CONVERTING LOCAL COORDINATES INTO GENOME-WIDE COORDINATES");
&convert_coordinates("genome.chunks.gff", "$output.cegma.local.gff", "$output.cegma.gff");


##########################################################################
# Produce summary report comparing predictions to set of 248 CEGs
##########################################################################
&print_header("EVALUATING RESULTS AND COMPARING TO SET OF 248 HIGHLY CONSERVED CEGS");
if ($cutoff_file eq "$CEGMA/data/profiles_cutoff.tbl") {
	my $options = "";
	$options = "-v " if ($verbose_flg);
	
    my $command = "completeness $options local_self.hmm_select.aln $CEGMA/data/completeness_cutoff.tbl > $output.completeness_report";
    &print_command("$command");
	system("$BIN/$command") && die "Can't run $command\n";
}


##########################################################################
##                                                                      ##
##                              ENDING                                  ##
##                                                                      ##
##########################################################################

## Removing extended files;
&clean_ext() unless $ext_flg;

## Removing temporary files
&clean_tmp() unless $temp_flg;

&print_header("CEGMA HAS FINISHED");



exit(0);








##########################################################################
##                                                                      ##
##                                                                      ##
##                                                                      ##
##                            SUBROUTINES                               ##
##                                                                      ##
##                                                                      ##
##                                                                      ##
##########################################################################

# Checking fasta format sequence files
sub check_fasta_format() {
    my $file = $_[0];
    my ($n, $c) = (undef, 0);
    open(TMP,"< $file");
    while (<TMP>) {
        next unless /^>/;
        />(\S+)\b/ && do {
            $n = $1;
            $c++;
            next unless $c>1;
        };
        &go_to_die("FATAL ERROR !!! Multiple locus names found.\n  File \'$file\' must contain only one sequence definition.");
    }
    &go_to_die("FATAL ERROR !!! There is no '>' line, locus name not found.\n  Please, verify your fasta file \'$file\'") unless defined($n);
    return $n;
}

# Deleting extended files
sub clean_ext() {
    # Obtaining the list of extended files
    my @files = qw( local.extended.tbl local_self.extended.tbl genome.blast genome.blast.gff 
       genome.chunks.fa genome.chunks.gff 
       local.geneid local.geneid.fa local.geneid.gff 
       local.genewise local.genewise.gff 
       local_self.geneid local_self.geneid.fa local_self.geneid.gff local_self.hmm_select.aln
		local_self.geneid.selected.gff
       local.geneid.selected.dna local.geneid.selected.fa 
       local.geneid.selected.gff local.geneid.selected.id local.hmm_select.aln 
       geneid_params/acc.logs geneid_params/intron.max geneid_params/self.param
       geneid_params/start.logs geneid_params/don.logs 
       geneid_params/coding.initial.5.logs 
       geneid_params/coding.transition.5.logs  );

   # Unlinking the temporary files if they exist
	Cegma::remove_files(@files);
    rmdir "geneid_params" if (-e "geneid_params");

}

# Deleting temporary files on TMP
sub clean_tmp() {
    # Obtaining the list of temporary files
    opendir(DIR,"$TMP");
    my @files = map { "$TMP/$_" } grep { /^$TMPROOT/ } readdir(DIR);
    closedir(DIR);
	Cegma::remove_files(@files);
}

# Checking input sequence files
sub exists_file() {
    my @files = @_;
    my ($n, $r) = (' ', 0);
    foreach $n (@files) {
        $r++ if (-e "$n");
    };
    return $r;
}

# Get a fixed length string from a given string and filling char/s.
sub fill_mid() { 
   my $l = length($_[0]); 
   my $k = int(($_[1] - $l)/2); 
   return ($_[2] x $k).$_[0].($_[2] x ($_[1] - ($l+$k)));
}

# writing die messages to STDERR and clean_tmp before exit.
sub go_to_die() { 
	(warn "@_\n") && &clean_tmp();
	&clean_ext() unless $ext_flg;
	print "\nEnding CEGMA\n\n";
	exit(1); 
}

# section headers to STDERR
sub print_header() { 
	return if ($quiet_flg);
	warn "\n\n";
	warn (("*" x 80)."\n** ".&fill_mid("@_",74," ")." **\n".("*" x 80)."\n");
	warn "\n";
}

# print the command to be executed
sub print_command() { 
	warn ("RUNNING: ".&fill_mid("@_",74," ")."\n") unless $quiet_flg;
}

# Run genome_map script
# this will run TBLASTN for proteins in each KOG and work out where the best 'chunks' in the genome are
# chunk = span in the genome where proteins in each KOG match

sub run_genome_map { 
    my ($genome_file, $prot_file, $tblastn_gff_file, $blastdb, $annotation, $interlen_val, $boundaries_val) = @_;

	my $file_prefix = "genome";
	
	# do we actually need to go any further?
	if (-s "$file_prefix.chunks.fa"){
		my $chunk_count = Cegma::count_seqs_in_fasta_file("$file_prefix.chunks.fa");
		warn "NOTE: $file_prefix.chunks.fa file already exists and contains $chunk_count candidate regions. Will use this file instead and skip running genome_map.\n" unless $quiet_flg;
		return;
	}
	
	# first check that FASTA headers don't include any examples with just digits before a space character
	my $digit_only_ids = `grep -cE \">[0-9]+( +|\$)\" $genome_file`;
	chomp($digit_only_ids);
	if($digit_only_ids){
		warn "There are $digit_only_ids sequences with FASTA headers that either contain only digits or have just digits followed by a space. E.g.\n";
		warn ">1\n";
		warn ">22 |\n";
		warn ">333 xyz\n\n";
		warn "These identifiers will break NCBI's blastdbcmd program and it won't be able to extract sequences later on from the BLAST database.\n";
		warn "Please reformat your FASTA headers.\n";		
		exit(1);
	}
	
	
	# set up options if some parameters have been defined
    my $options = "-n $file_prefix -p $prot_num -o $interlen_val -c $boundaries_val -t $threads ";
    ($options .= " -b $tblastn_gff_file ") if ($tblastn_gff_file);	
    ($options .= " -d $blastdb ")          if ($blastdb);	
    ($options .= " -a $annotation ")       if ($annotation);
	($options .= " -v ")                   if ($verbose_flg);
    
    my $genome_map_command = "genome_map  $options $prot_file $genome_file 2>$output.cegma.errors";

    &print_command("$genome_map_command");
	system("$BIN/$genome_map_command") && &go_to_die("FATAL ERROR when running genome_map $?: \"$!\"\n");

	# did we end up with a non-zero byte chunks file?
	prematurely_end_program() if (-z "$file_prefix.chunks.fa");
}

sub run_local_map { 
    my ($outputfile, $option) = @_;
	my  $options = "-n $outputfile $option -h $hmm_profiles -i $hmm_prefix ";
	$options .= "-v " if ($verbose_flg);
    my $local_map_command = "local_map $options genome.chunks.fa 2>$output.cegma.errors";

    &print_command("$local_map_command");

    system("$BIN/$local_map_command") && &go_to_die("FATAL ERROR when running local map $?: \"$!\"");

	if (-z "local.genewise.gff"){
	   &go_to_die("\nGenewise did not produce any suitable alignments. I.e. no orthologous genes were found in the regions identified by TBLASTN. Aborting CEGMA.\n\n");
	}

	my $geneid_count = Cegma::count_seqs_in_fasta_file("$outputfile.geneid.fa");
	warn "NOTE: created $geneid_count geneid predictions\n" unless ($quiet_flg);	
	prematurely_end_program() if ($geneid_count == 0);
}

 
 sub prematurely_end_program{
	&print_header("EXITING CEGMA PREMATURELY");

	print "CEGMA is stopping because there are no core gene remaining in the data at this stage of analysis.\n\n\n";
	## Removing extended files;
	&clean_ext() unless $ext_flg;

	## Removing temporary files
	&clean_tmp() unless $temp_flg;
	
	exit(1);
 }

# convert the coordinates of the fragments to the genomic ones
sub convert_coordinates() {
    my ($file1, $file2, $file3) = @_;
    my (%pos1, %pos2, %strand, %seq);
    
    # first we read the genomic fragments (chunks)
    open(FILE1, "<", "$file1") or die "Can't read from $file1\n";
    while (<FILE1>) {
		my @gff = split;
		$seq{$gff[8]}    = $gff[0];
		$pos1{$gff[8]}   = $gff[3];
		$pos2{$gff[8]}   = $gff[4];
		$strand{$gff[8]} = $gff[6];
    }
    close(FILE1);

    open(FILE2, "<", "$file2") or die "Can't read from $file2\n";;
    open(FILE3, ">", "$file3") or die "Can't write to $file3\n";;
    while (<FILE2>) {
		next if (m/^#/);
		my @gff = split;

		if ($strand{$gff[0]} eq "+"){
	    	printf FILE3  "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $seq{$gff[0]}, $PROGRAM, $gff[2], 
	            $gff[3] + $pos1{$gff[0]} - 1, $gff[4] + $pos1{$gff[0]} - 1, $gff[5], "+", $gff[7], $gff[0];
		} elsif ($strand{$gff[0]} eq "-"){
	    	printf FILE3  "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $seq{$gff[0]}, $PROGRAM, $gff[2], 
	            $pos2{$gff[0]} - $gff[4] + 1, $pos2{$gff[0]} - $gff[3] + 1, $gff[5], "-", $gff[7], $gff[0];
		} 
    }

    close(FILE2);
    close(FILE3);

    return(0);
}

# get a fixed length string from a given string and filling char/s.
sub Which_Options() {
    GetOptions( 
             "g|genome=s"      => \$genome       , # genome fasta file
             "p|protein=s"     => \$protein      , # protein fasta file
	         "o|output=s"      => \$output       , # output file name
             "c=f"             => \$score_cutoff , # tblastx score cutoff
             "s|genewise=s"    => \$genewise_file, # shrink hsp's by
             "t|tblastn_gff=s" => \$tblastn_gff  , # read tblastn results from existing GFF file
             "d|blastdb=s"     => \$blastdb      , # blastdb file
	         "a|annotation=s"  => \$annotation   , # annotation file
             "vrt"             => \$vrt_flg      , # vertebrate
             "mam"             => \$mam_flg      , # mamalian
             "hmm_profiles=s"  => \$hmm_profiles , # directory to find the hmms
             "hmm_prefix=s"    => \$hmm_prefix   , # prefix for the profiles
             "prot_num=f"      => \$prot_num     , # number of protiens per profile
             "cutoff_file=s"   => \$cutoff_file  , # file with the cutoffs
	         "max_intron=f"    => \$max_intron   , # max intron for geneid
             "v|verbose"       => \$verbose_flg  , # verbose    
             "q|quiet"         => \$quiet_flg    , # quiet mode
             "ext"             => \$ext_flg      , # kepp extended files
             "tmp"             => \$temp_flg     , # keep temporary files
             "h|help|?"        => \$help_flg     , # print help
             "T|threads=i"     => \$threads      , # number of threads
                ) || &print_help();
    &print_help if $help_flg;

   # cheking if files exists
   &go_to_die("FATAL ERROR!!! --genome or --blastdb must be defined.") 
	   unless (defined($genome) || defined($blastdb));
   
   my $file_number = &exists_file($genome);
   &go_to_die("FATAL ERROR! Input file: \"$genome\" does not exist.")
	   unless $file_number == 1; 
     

   &go_to_die("FATAL ERROR!!! --verbose and --quiet are incompatible options.")
	   if ($verbose_flg && $quiet_flg);

   &go_to_die("FATAL ERROR!!! --intron_max should be bigger than 50 bp.")
	   if ($max_intron && $max_intron<50);
   
   # default assigments 
   $protein = "$CEGMA/data/kogs.fa" if (!(defined($protein)));
   $output = "output" if (!(defined($output)));
   $tblastn_gff = 0   if (!(defined($tblastn_gff)));
   $annotation = 0   if (!(defined($annotation)));
   $hmm_profiles = "$CEGMA/data/hmm_profiles"   if (!(defined($hmm_profiles)));
   $cutoff_file = "$CEGMA/data/profiles_cutoff.tbl" if (!(defined($cutoff_file))); 
   $prot_num = 6 if (!(defined($prot_num)));
   $hmm_prefix = "KOG" if (!(defined($hmm_prefix)));
   #  $GeneidParam   = "$geneid_param" if defined($geneid_param);
   $verbose_str = " -v " if $verbose_flg;
   $interlen = "20000" if $vrt_flg;
   $boundaries = "5000" if $vrt_flg;
   $interlen = "40000" if $mam_flg;
   $boundaries = "5000" if $mam_flg;
   $threads = 2 if (!defined($threads));
   
} 

# prints help 
sub print_help() {
    open(HELP, "| cat") ;
    print HELP <<"+++EndOfHelp+++";

                              $PROGRAM


PROGRAM:
                        $PROGRAM - $VERSION

                Core Eukaryotic Genes Mapping Approach

USAGE:
        
    $PROGRAM [options] <-g genomic_fasta_sequence> 

DESCRIPTION:

    CEGMA (Core Eukaryotic Genes Mapping Approach) is a pipeline for
    building a set of high reliable set of gene annotations in
    virtually any eukaryotic genome. It combines tblastn, genewise,
    hmmer, with geneid, an "ab initio" gene prediction program.

REQUIRES:

    $PROGRAM needs the pre-installation of the following software:

     - geneid (geneid v 1.4)
     http://genome.imim.es/software/geneid/
     - genewise (wise2.2.3-rc7, which also requires glib to be installed)
     http://www.ebi.ac.uk/Wise2/
     - hmmer (HMMER 3.0) 
     http://hmmer.janelia.org/
     - NCBI blast+ (TBLASTN 2.2.25 [31-Mar-2011])
     ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
 
    Check that you have the right version of the previous software.
    If you are using a different version and experience any problems,
    please let us know !!!

ENVIRONMENT VARIABLES:

    You can specify the path where $PROGRAM can find the default files
    with the shell variable \"$PROGRAM\". 

    You also can specify the path for the tempotary files with the
    shell variable \"CEGMATMP\". Default value is /tmp.

	Setting those vars in Bourne-shell and C-shell:

     o Using a Bourne-Shell (e.g. bash):
           export CEGMA="path"
           export CEGMATMP="path"
           export PERL5LIB="\$PERL5LIB:\$CEGMA/lib" 

     o Using a C-Shell:
           setenv CEGMA "path"
           setenv CEGMATMP "path"
           setenv PERL5LIB "\$PERL5LIB:\$CEGMA/lib"

 	Genewise will also require that you set the \$WISECONFIGDIR environment variable.
   
 
COMMAND-LINE OPTIONS:

	 Available options and a short description are listed here:

     -g, --genome      fasta file of the query sequence.

     -p, --protein     fasta file of the protein sequences.
                        (only used to run a subset of the 458 CEGs)

     -o, --output      ouput file prefix.

     -d, --blastdb     blast database for the genome sequence.

     -t, --tblastn_gff  gff file containing the tblastn results
                        (this option skips the blast step).

     -s, --genewise    gff file with the genewise alignment coordinates
                        (this options skips the blast step).

     -a, --annotation  gff file with gene annotations
                        (computes the coordinates for the annotations).

     --vrt             Optimization for vertebrate genomes
                        (intron size 20,000 bp).

     --mam             Optimization for mamalian genomes
                        (intron size 40,000 bp).

     --max_intron      Max intron size.

     -T, --threads     Specify number of processor threads to use
	 
     -v, --verbose     Verbose mode, shows progress of each KOG as it is processed

     -q, --quiet       Quiet mode, do not show any message/warning                       

     -h, --help        Shows this help.

     --ext             Extended output which keeps all intermediate files.

     --tmp             Keep temporary files.

     These variables can be used to run other than default HMM profiles:
     
     --prot_num        Number of proteins in the fasta file 
                       They have to be in consecutive order in the fasta file
                          (default: 6)
     --cutoff_file     File with the cutoff for the HMMER alignments
                          (default: \$CEGMA/data/profiles_cutoff.tbl) 
     --hmm_prefix      Each protein ID must have "___" followed by the hmmprefix 
                       and a number (ex. At3g02190___KOG1762)
                           (default: KOG)
     --hmm_directory   Directory that contains the hmm files. The files must be
                       named hmm_prefix(number).hmm  ex. KOG1762.hmm
                          (default: \$CEGMA/data/hmm_profiles)     

BUGS:    
  
    Please report any problem to 'korflab\@ucdavis.edu'.

AUTHOR:  

    $PROGRAM has been developed by Genis Parra, Keith Bradnam, and Ian Korf.



GNU-GPL (C)                      May 2014                         $PROGRAM
+++EndOfHelp+++
    close(HELP);
    exit(1);
} # print_help

