#!/usr/bin/perl 

##########################################################################
#                                                                        #
#                              hmm_select                                #
#                                                                        #
#                         CEGMA pipeline code                            #
#                                                                        #
##########################################################################

use strict; use warnings; use sigtrap;
use FAlite;
use Cegma;
use Getopt::Std;

use vars qw($opt_p $opt_o $opt_i $opt_v $opt_t);

getopts('p:o:i:vt:');

my $PREFIX = "output";
my $PROTNUM;
my $HMM_PREFIX ="KOG";
my $THREADS = 2;

$PREFIX     = $opt_o if $opt_o; # output prefix for filenames
$PROTNUM    = $opt_p if $opt_p;
$HMM_PREFIX = $opt_i if $opt_i;
$THREADS    = $opt_t if $opt_t;

die "

# This script takes a protein multifasta file and run it against a hmm 
# profile using hmmsearch. 

hmm_select.pl <hmm_directory> <protein_fasta_file> <cutoff_file>
options:
   -p <protein file>
   -i <hmm prefix, e.g. 'KOG'>
   -o <output file prefix>
   -t <number of threads>
   -v verbose mode, show progress of each KOG

" unless (@ARGV == 3);

my ($hmm_directory, $prot_file, $cutoff_file) = @ARGV;

# Security check to not overwrite the results
die "File does not exist: $prot_file \n" if (not -e $prot_file);

# look up pre-generated score cutoffs for each KOG
my %score;
my %score_cutoff;
my %best_pred;

open(CUTOFF_FILE, "<", "$cutoff_file") or die "Can't read from $cutoff_file\n";

while (<CUTOFF_FILE>){
    if (/($HMM_PREFIX[0-9]+)\s(.*)/){
		$score_cutoff{$1} = $2;
		$score{$1} = 0;
    }
}
close(CUTOFF_FILE);

my $protein_count = Cegma::count_seqs_in_fasta_file($prot_file);
my $protein_counter = 0;

open(EXTENDED_FILE, ">", "$PREFIX.extended.tbl") or die "Can't write to $PREFIX.extended.tbl\n";

# Opening the protein fasta files
open(PROT_FILE, "<", $prot_file) or die "Can't read from $prot_file\n";
my $PROT = new FAlite(\*PROT_FILE);

# remove the hmm_select output file if present

if (-e "$PREFIX.hmm_select.aln"){
	warn "WARNING: $PREFIX.hmm_select.aln already exists. Removing\n";
	unlink ("$PREFIX.hmm_select.aln");
}


# Main loop of execution
# looping over protein sequences of genes predicted by geneid
# can still have more than one protein per KOG at this stage

while (my $prot_entry = $PROT->nextEntry) {

    # Opening temporary fasta files
    open(TMP_PROT, ">", "/tmp/prot$$.fa") or die "Can't write to /tmp/prot$$.fa\n";
    Cegma::print_fasta($prot_entry,\*TMP_PROT);
    my ($kog_id)      = $prot_entry->def =~ /($HMM_PREFIX\d+)/;
    my ($full_kog_id) = $prot_entry->def =~ /($HMM_PREFIX\d+\.\d+)/;

    close(TMP_PROT);

	$protein_counter++;    
    print "Processing geneid prediction $protein_counter/$protein_count: $full_kog_id\n" if ($opt_v);

    # running hmmsearch
    system ("hmmsearch --cpu $THREADS $hmm_directory/$kog_id.hmm /tmp/prot$$.fa > /tmp/hmm$$.output") && die "Can't run hmmsearch\n";

	# want to find geneID predictions that are above threshold score
    open(TMP_HMM, "<", "/tmp/hmm$$.output") or die "Can't read from /tmp/hmm$$.output\n";
    while (<TMP_HMM>){
		next if m/^#/;
		
		# match score from hmmsearch output
		if (/\s+\S+\s+([\d\.]+).*$full_kog_id/){
			my $score = $1;
	    	print EXTENDED_FILE "$full_kog_id $score $score_cutoff{$kog_id}\n";
			# only want the best score when there are multiple predictions per KOG
	    	if ($score > $score_cutoff{$kog_id} && $score > $score{$kog_id}){
				$score{$kog_id} = $score;
				$best_pred{$kog_id} = $full_kog_id;
	    	}
			last;
		}
    }
    close(TMP_HMM);

    # copy the output in the final file:
    system ("cat /tmp/hmm$$.output >> $PREFIX.hmm_select.aln") && die "Can't concatenate hmm$$.output to $PREFIX.hmm_select.aln\n";

}
close (PROT_FILE);
close (EXTENDED_FILE);



# Print the geneid proteins that pass the threshold scores, just one protein for each KOG 
# (the highest scoring when multiple predictions exceed threshold). Also want to print IDs to a file
# and find the matching ID from the main GFF file to add to a selected GFF file, and then
# grab the corresponding sequence from genome.chunks.fa

if (-s "$PREFIX.geneid.selected.gff"){
	warn "WARNING: $PREFIX.geneid.selected.gff already exists. Removing\n";
	unlink("$PREFIX.geneid.selected.gff");	
}

open(PROT_FILE, "<", $prot_file) or die "Can't read from $prot_file\n";
my $PROT2 = new FAlite(\*PROT_FILE);

open(my $id_out,    ">", "$PREFIX.geneid.selected.id") or die "Can't write to $PREFIX.geneid.selected.id\n";
open(my $fasta_out, ">", "$PREFIX.geneid.selected.fa") or die "Can't write to $PREFIX.geneid.selected.fa\n";

# keep track of the selected ID
my %good_ids;

while (my $prot_entry = $PROT2->nextEntry) {
    my ($kog_pred_id) = $prot_entry->def =~ /($HMM_PREFIX\d+.[0-9]+)/;
    my ($kog_id) = $prot_entry->def =~ /($HMM_PREFIX\d+)/;
    if ($best_pred{$kog_id} && $best_pred{$kog_id} eq $kog_pred_id){
		Cegma::print_fasta ($prot_entry, $fasta_out);
		print $id_out "$kog_pred_id\n";
		$good_ids{$kog_pred_id} = 1;
		system ("grep -w $kog_pred_id $PREFIX.geneid.gff >> $PREFIX.geneid.selected.gff");
    }
}
close($id_out);
close($fasta_out);



# can now look up good IDs to extract sequence from genome.chunks.fa 
my $infile  = "genome.chunks.fa";
my $outfile = "$PREFIX.geneid.selected.dna"; 

open(my $in,  "<", $infile)  or die "Can't read from $infile\n";
open(my $out, ">", $outfile) or die "Can't write to $outfile\n";

my $genome = new FAlite(\*$in);
while (my $entry = $genome->nextEntry) {
	my $header = $entry->def;
	my ($id) = $header =~ m/>(\S+)\s+/;
	if(exists $good_ids{$id}){
		# just need the KOG ID in the header for this file
		$entry =~ s/(>\S+)\s+.*/$1/;
		Cegma::print_fasta($entry, $out);
	} 
}
close($in);
close($out);


exit;

