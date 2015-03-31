#!/usr/bin/perl

##########################################################################
#                                                                        #
#                              genome_map                                #
#                                                                        #
#                         CEGMA pipeline code                            #
#                                                                        #
##########################################################################

use strict; use warnings; use sigtrap;
use FAlite;
use Cegma;
use Getopt::Std;
use vars qw($opt_p $opt_n $opt_b $opt_d $opt_a $opt_o $opt_c $opt_v $opt_t);

##########################################################################
##                                                                      ##
##                     VARIABLES  AND DEFINITIONS                       ##
##                                                                      ##
##########################################################################
 
my $CEGMA = defined($ENV{'CEGMA'}) ? $ENV{'CEGMA'} : '.';
my $TMP   = defined($ENV{'CEGMATMP'}) ? $ENV{'CEGMATMP'} : '/tmp';

getopts('p:n:b:d:a:o:c:vt:');

my $PREFIX = "";
my $PROTNUM;
my $TBLASTN_GFF_FILE;
my $ANNOT_FILE;
my $BLASTDB = "$TMP/genome$$.blastdb";
my $OVERLAP;
my $BOUNDARIES;
my $THREADS = 2;

$PREFIX           = $opt_n if $opt_n;
$PROTNUM          = $opt_p if $opt_p;
$BLASTDB          = $opt_d if $opt_d;
$TBLASTN_GFF_FILE = $opt_b if $opt_b;
$ANNOT_FILE       = $opt_a if $opt_a;
$OVERLAP          = $opt_o if $opt_o;
$BOUNDARIES       = $opt_c if $opt_c;
$THREADS          = $opt_t if $opt_t;
die "

# This script takes a protein multifasta file and a genomic multifasta file
# and runs blast to get the best n candidates.

# If a blast file is provided it gets the data from there.

# we considered the sequence of the proteins from the same cluster to be 
# in  order.

genome_map.pl <protein_fasta_file> <genomic_fasta_file>
   -n <prefix of the output files>
   -p <number of species proteins per local proteins>
   -o overlap : allowable distance between separate HSPs that are to be considered 1 'chunk'
   -c boundaries : distance subtracted/added to final coordinates of each chunk
options:
   -b <tblastn gff file> 
   -d <genomic blast database>
   -a <annotation_gff>
   -v verbose mode - prints information for each KOG that is processed
   -t <number of threads>

" unless (@ARGV == 2 && $opt_p && $opt_o && $opt_c);

my ($prot_file, $genome_file) = @ARGV; 

##########################################################################
##                                                                      ##
##                          INITIAL CHECKINGS                           ##
##                                                                      ##
##########################################################################

# check command-line options (only checking for one file as part of blast database)
die "genome_map: File does not exist: $prot_file \n"                      if (not -e $prot_file);
die "genome_map: File does not exist: $genome_file \n"                    if (not -e $genome_file);
die "genome_map: The specified BLAST database does not exist: $opt_d\n"   if (defined $opt_d && not -e "$opt_d.nsq");
die "genome_map: The specified TBLASTN GFF file does not exist: $opt_b\n" if (defined $opt_b && not -e "$opt_b");


#if (!(-e "$ANNOT_FILE" )) {
#    print "File does not exist: $ANNOT_FILE \n";
#    exit(1);
#}

# do we have an existing files that will be overwritten?
print "WARNING: $PREFIX.blast already exists, this will be overwritten\n"      if (-e "$PREFIX.blast");
print "WARNING: $PREFIX.blast.gff already exists, this will be overwritten\n"  if (-e "$PREFIX.blast.gff");
print "WARNING: $PREFIX.chunks.gff already exists, this will be overwritten\n" if (-e "$PREFIX.chunks.gff");


##########################################################################
##                                                                      ##
##                          COMPUTING DATA                              ##
##                                                                      ##
##########################################################################

# make BLAST database (unless a pre-existing database has already been specified)
if (not defined($opt_d)){
	my $command = "makeblastdb -in $genome_file -dbtype nucl -parse_seqids -out $BLASTDB";
	system($command) && die "Can't run $command\n";	
} else{
	print "NOTE: Using specified BLAST database: $opt_d\n";
}

print "NOTE: Using specified TBLASTN GFF file ($TBLASTN_GFF_FILE) instead of running BLAST\n" if (defined $TBLASTN_GFF_FILE);

# Opening the protein fasta files
open(PROT_FILE, "<", $prot_file) or die "Can't read from $prot_file\n";
my $PROT = new FAlite(\*PROT_FILE);

# Main loop of execution
# Need to write all six proteins from each KOG to separate FASTA file and run TBLASTN search
# with just those sequences. Then convert output into GFF file

my $kog_counter = 0;

while (my $prot_entry = $PROT->nextEntry) {
	
    # Opening temporary fasta files
    open(TMP_PROT, "> $TMP/prot$$.fa") or die "Can't write to $TMP/prot$$.fa\n";
    Cegma::print_fasta($prot_entry,\*TMP_PROT);
    my ($kog_id) = $prot_entry->def =~ /^>\S+___(\S+)/;

    for (my $c = 1; $c < $PROTNUM; $c++){
		$prot_entry = $PROT->nextEntry;
		Cegma::print_fasta($prot_entry,\*TMP_PROT);	
    }
    close(TMP_PROT);

	$kog_counter++;
	print "Processing KOG #$kog_counter: $kog_id \n" if ($opt_v);

    # no need to run TBLASTN if a GFF output file already exists
    if (!($TBLASTN_GFF_FILE)) {
		# running blast
		system ("tblastn -db $BLASTDB -query $TMP/prot$$.fa \\
	             -word_size 6 -max_target_seqs 5 -seg yes -num_threads $THREADS -lcase_masking \\
				 -outfmt \"7 sseqid sstart send sframe bitscore qseqid\"\\
	             > $TMP/blast$$") && die "Can't run tblastn\n";

		system ("cat $TMP/blast$$ >> $PREFIX.blast") && die "Can't concatenate $TMP/blast$$ to $PREFIX.blast\n";

		# now create a GFF file from the BLAST output
		open(BLAST, "<", "$TMP/blast$$")     or die "Can't read from $TMP/blast$$";
		open(GFF,   ">", "$TMP/blast$$.gff") or die "Can't create blast$$.gff";

		while(<BLAST>){
			chomp;
			next if m/^#/; 
			my ($s, $s_start, $s_end, $sframe, $score, $q) = split;			

			# tidy up subject frame (if on negative strand)
			$sframe =~ s/^-//;
			
			# determine strand and print appropriate information
			if ($s_start > $s_end){
				print GFF "$s\tTBLASTN\thsp\t$s_end\t$s_start\t$score\t-\t$sframe\t$q\n";				
			} else{
				print GFF "$s\tTBLASTN\thsp\t$s_start\t$s_end\t$score\t+\t$sframe\t$q\n";
			}

		}
		close(BLAST);
		close(GFF);

		# want to have GFF file sorted by 1) seq ID and 2) start coord
		system("sort -nk 1 -nk 4 $TMP/blast$$.gff > $TMP/blast$$.gff2") && die "Can't sort $TMP/blast$$.gff";
		rename("$TMP/blast$$.gff2", "$TMP/blast$$.gff")             or die "Can't rename sorted GFF file\n";

		system ("cat $TMP/blast$$.gff >> $PREFIX.blast.gff") && die "Can't concanenate $TMP/blast$$.gff to $PREFIX.blast.gff\n";
    } else {
		# if existing TBLASTN GFF file is specified can just grab information from that
		system ("grep $kog_id $TBLASTN_GFF_FILE > $TMP/blast$$.gff") && die "Can't grep KOG IDs from $TBLASTN_GFF_FILE\n";
		system("sort -nk 1 -nk 4 $TMP/blast$$.gff > $TMP/blast$$.gff2") && die "Can't sort $TMP/blast$$.gff";
		rename("$TMP/blast$$.gff2", "$TMP/blast$$.gff")             or die "Can't rename sorted GFF file\n";

    }
	
	# convert GFF file for current KOG to a 'chunks' GFF file which shows the top 5 scoring regions
	# of the input sequence
	convert_gff_to_chunks("$TMP/blast$$.gff", $kog_id);
#	exit;
}

close (PROT_FILE);

my $chunk_count = Cegma::count_seqs_in_fasta_file("$PREFIX.chunks.fa");
warn "Found $chunk_count candidate regions in $PREFIX.chunks.fa\n" if ($opt_v);


# removing temporary files
Cegma::remove_files("$TMP/prot$$.fa", "$TMP/chunks$$.gff", "$TMP/blast$$", "$TMP/blast$$.gff", "$TMP/chunks$$.fa", glob("$TMP/genome$$.blastdb.n*"));


exit(0);

##########################################################################
##                                                                      ##
##                              THE END                                 ##
##                                                                      ##
##########################################################################





##########################################################################
##                                                                      ##
##                            SUBROUTINES                               ##
##                                                                      ##
##########################################################################




# takes a GFF file which contains TBLASTN output of all proteins in 1 KOG searched
# against input sequence. Find the top 5 scoring regions (individual HSPs are allowed
# to overlap by $OVERLAP) in the file and outputs these to a new GFF 'chunks' file.
# coordinates of each chunk in the chunks file are extended by $BOUNDARIES

sub convert_gff_to_chunks{
	
	my ($gff_file, $kog_id) = @_;
	
	
	# for each GFF file we want to track the start and end of each chunk that we find
	# as well as keep track of the IDs, scores and total counts of all chunks
	my ($line_count, $chunk_count) = (0, 0);
	my ($chunk_seq, $chunk_min, $chunk_max, $chunk_strand, $chunk_frame, $chunk_score); 

	# track whether current line of GFF file belongs to same chunk or not
	my $same_chunk = 0; 
	
	# hash to track chunks, stores all basic chunk information tied to chunk ID number as key
	# six secondary hash keys are: min, max, score, frame, strand, seq
	my %chunks;
	open(GFF, "<", "$gff_file") or die "can't open $gff_file";

	while(<GFF>){
		$line_count++;
		my ($seq, $source, $feature, $start, $end, $score, $strand, $frame, undef) = split;			

		# special case for first line
 		if ($line_count == 1){
			$chunk_count++;
			($chunk_seq, $chunk_min, $chunk_max, $chunk_strand, $chunk_frame, $chunk_score) = ($seq, $start, $end, $strand, $frame, $score);
			$chunks{$chunk_count}{min}    = $start;
			$chunks{$chunk_count}{max}    = $end;
			$chunks{$chunk_count}{score}  = $score;
			$chunks{$chunk_count}{strand} = $strand;
			$chunks{$chunk_count}{seq}    = $seq;
			$chunks{$chunk_count}{frame}  = $frame;			
			next;
		}


		# looking for overlapping matches to same sequence on same strand
		if (($seq eq $chunk_seq) && ($strand eq $chunk_strand)){
			# same coords
			if (($start == $chunk_min) && ($end == $chunk_max)){
				($chunk_score += $score);
				$same_chunk = 1;				
			} 
			# left overhang
			elsif (($start >= ($chunk_min - $OVERLAP)) && ($start <= $chunk_min) && ($end >= ($chunk_max - $OVERLAP)) && ($end <= $chunk_max)){
				$chunk_score += $score;
				$chunk_min = $start;
				$same_chunk = 1;
			}
			# right overhang
			elsif (($start <= ($chunk_min + $OVERLAP)) && ($start >= $chunk_min) && ($end <= ($chunk_max + $OVERLAP)) && ($end >= $chunk_max)){
				$chunk_score += $score;
				$chunk_max = $end;
				$same_chunk = 1;
			} 
			# double overhang
			elsif (($start >= ($chunk_min - $OVERLAP)) && ($start <= $chunk_min) && ($end <= ($chunk_max + $OVERLAP)) && ($end >= $chunk_max)){
				$chunk_score += $score;
				$chunk_min = $start;
				$chunk_max = $end;
				$same_chunk = 1;				
			} 
			# double underhang
			elsif (($start <= ($chunk_min + $OVERLAP)) && ($start >= $chunk_min) && ($end >= ($chunk_max - $OVERLAP)) && ($end <= $chunk_max)){
				$chunk_score += $score;
				$same_chunk = 1;				
			} 
			# if here then we are on the same strand of the same sequence, but don't overlap
			else{
				$same_chunk = 0;
			}					
		} 
		# if we get here then we are not on the same sequence and/or strand
		else{
			$same_chunk = 0;
		}
			

		if ($same_chunk == 0){
			($chunk_seq, $chunk_min, $chunk_max, $chunk_strand, $chunk_score) = ($seq, $start, $end, $strand, $score);
			$chunk_count++;
		} 
		
		# first place where we can update %chunks
		$chunks{$chunk_count}{min}    = $chunk_min;
		$chunks{$chunk_count}{max}    = $chunk_max;
		$chunks{$chunk_count}{score}  = $chunk_score;
		$chunks{$chunk_count}{strand} = $chunk_strand;
		$chunks{$chunk_count}{seq}    = $chunk_seq;
		$chunks{$chunk_count}{frame}  = $chunk_frame;
	}

	close(GFF);


	# want to show just the top 5 chunks and write these to a new file
	open (CHUNKS, ">", "$TMP/chunks$$.gff");
	
	my $counter = 0;
	foreach my $id (sort {$chunks{$b}{score} <=> $chunks{$a}{score}} keys %chunks){
		last if $counter == 5;

		# modify start and end of chunk by $BOUNDARIES
		my $start      = $chunks{$id}{min} - $BOUNDARIES;		
		my $end        = $chunks{$id}{max} + $BOUNDARIES;
		my $strand     = $chunks{$id}{strand};
		my $seq        = $chunks{$id}{seq};
		my $new_kog_id = "$kog_id.$id";

		# just make sure we don't end up with negative coordinates
		$start = 1 if $start < 1;

		# NOT SURE WHAT THIS PART IS DOING???
		# Getting the genomic coordinates and converted to local
		# gawk script to get a reverse the coordinates
		if (defined ($ANNOT_FILE)){
		    system ("awk '\$1 == \"$seq\" && \$4 > $start && \$5 < $end' \\
                      $ANNOT_FILE \\
                    | $CEGMA/bin/local_coord $start $end $strand $new_kog_id \\
                    >> $PREFIX.chunks.ann.gff") && die "Can't run awk command\n";
		}


		# now want to check that expanded $end coordinate might not exceed length of available sequence
		# form command to extract sequence from BLASTDB
		
		# first want length of target sequence
		my $blastdbcmd_command = "blastdbcmd -dbtype nucl -outfmt %l -db $BLASTDB -entry $seq";
		my $seq_length = `$blastdbcmd_command`;
		chomp($seq_length);
		
		# if length is greater than $end, modify $end
		$end = $seq_length if ($end > $seq_length);
		
		# now extract sequence
		$blastdbcmd_command = "blastdbcmd -dbtype nucl -range $start-$end -db $BLASTDB -entry $seq";
		$blastdbcmd_command .= " -strand minus" if $strand eq "-";
#		print "running: $blastdbcmd_command\n";
		system("$blastdbcmd_command | sed 's/>\\(.*\\)/>$new_kog_id \\1/' > $TMP/chunks$$.fa") && die "Can't run $blastdbcmd_command\n";


		# as we add some context to generate the candidate region, these coordinates
	   	# can exceed the length of genomic fragments and that can be a problem when
	   	# converting the final predictions to genomic coordinates again (on the negative strand)
	   	# we are going to check the length of the fragments and correct the coordinates 
	   	# to build the final gff file. 
	    open(CHUNKS_SEQ, "<", "$TMP/chunks$$.fa") or die "Can't read from chunks$$.fa\n";
		my $CHUNK = new FAlite(\*CHUNKS_SEQ);
		my $chunk_entry = $CHUNK->nextEntry;
	    my $chunk_length = length($chunk_entry->seq);
	    if ($end - $start + 1 != $chunk_length){
#			print "In the place where bad things happen\n";
	    	$end = $start + $chunk_length - 1;
	    }
		close(CHUNKS_SEQ);
		
		# can now print final chunk result to file
		print CHUNKS "$seq\tTBLASTN\thsp\t$start\t$end\t$chunks{$id}{score}\t$strand\t$chunks{$id}{frame}\t$new_kog_id\t\n";

	    system ("cat $TMP/chunks$$.fa >> $PREFIX.chunks.fa;") && die "Can't concatenate $TMP/chunks$$.fa to $PREFIX.chunks.fa\n";
	
		$counter++;
	}
	close CHUNKS;
    system ("cat $TMP/chunks$$.gff >> $PREFIX.chunks.gff") && die "Can't concatenate $TMP/chunks$$.gff to $PREFIX.chunks.gff\n";

	return;
}
