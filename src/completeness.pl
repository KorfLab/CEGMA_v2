#!/usr/bin/perl

##########################################################################
#                                                                        #
#                             completness                                #
#                                                                        #
#                         CEGMA pipeline code                            #
#                                                                        #
##########################################################################

use strict; use warnings; use sigtrap;
use FAlite;
use Getopt::Std;
use List::Util qw(max min);
use vars qw($opt_e $opt_v);

getopts('ev');

die "

# This script takes a hmmsearch output file and a cutoff paramaters file 
# and give the statistic of completeness

completeness <hmmsearch_output> <cutoff_file>

options:
 -v           verbose (will print details of 248 CEGs missing from genome)
 -e           alternative length cutoffs (0%, 50%, 60%, 70%, 80%, and 90%)
" unless (@ARGV == 2);


# want to process just the one top scoring gene prediction in each KOG
my ($hmm_file, $cutoff_file) = @ARGV;


######################################################################
#
#     Security check
#
######################################################################

if (!(-e "$hmm_file" )) {
    print "File does not exist: $hmm_file \n";
    exit(1);
}


if (!(-e "$cutoff_file" )) {
    print "File does not exist: $cutoff_file \n";
    exit(1);
}


######################################################################
#
#      Reading the cutoff and conservation group for each protein
#
######################################################################

my %score;
my %conservation;
my %conservation_groups;
my %score_cutoff;
my %prot_length;
my $total_proteins;

my %length_cutoff;

if ($opt_e) {
   %length_cutoff = (
		     "Cutoff 90%" => '90',
		     "Cutoff 80%" => '80',
		     "Cutoff 70%" => '70',
		     "Cutoff 60%" => '60',
		     "Cutoff 50%" => '50',
		     "Cutoff 0%" => '0',
		     );

} else {
   %length_cutoff = (
		     Complete => '70',
		     Partial  => '0',
		     );
}

open(CUTOFF_FILE, "<", "$cutoff_file") or die "Can't read from $cutoff_file\n";
while (<CUTOFF_FILE>){
    if (/(.*) (.*) (.*) (.*)/){
		$conservation{$1} = $2;
		$conservation_groups{$2}++;
		$score_cutoff{$1} = $3;
		$prot_length{$1} = $4;
		$score{$1} = 0;
		$total_proteins++;
    }    
}

close(CUTOFF_FILE);


######################################################################
#
#     Reading the hmmsearch output and getting the score and length 
#     of the alignments
#
######################################################################

my %over_cutoff;
my %over_cutoff_locus;
my %over_cutoff_group;
my %over_cutoff_group_locus;
my $locus;
my $warning_version = 0;

open(HMMSEARCH, "<", "$hmm_file") or die "Can't read from $hmm_file\n";

# check HMMER version
<HMMSEARCH>;
my $second_line = <HMMSEARCH>;
my ($version) = $second_line =~ m/HMMER (\d+\.\d+)/;
$warning_version = 1 if ($version ne "3.0");

while (<HMMSEARCH>){

	next if m/^#/;		
	if (m/^Query:/){

		# first grab KOG
		# Query:       KOG0019  [M=709]
		my ($kog) = $_ =~ m/Query:\s+(\S+?)\s+\[/;
		
		# skip next 4 lines
		#Scores for complete sequences (score includes all domains):
		#   --- full sequence ---   --- best 1 domain ---    -#dom-
		#    E-value  score  bias    E-value  score  bias    exp  N  Sequence                                           Description
		#    ------- ------ -----    ------- ------ -----   ---- --  --------                                           -----------
		<HMMSEARCH>; <HMMSEARCH>; <HMMSEARCH>; <HMMSEARCH>; 

		# grab specific KOG ID and score from next line
		#   1.6e-118  382.3  21.0    3.3e-63  199.1   4.6    3.0  3  KOG0019.6_1|geneid_v1.4_predicted_protein_1|682_AA
		my $line = <HMMSEARCH>;
		my ($score, $kog_id) = $line =~ m/\s+[\d\.e-]+\s+([\d\.]+).*(KOG\d+\.\d+)_/;

		#		next unless exists($best_kog_ids{$kog_id});
		
		# some geneID predictions will have produce no significant HMMER
		# hits, and you won't be able to extract KOG ID, so just skip
		next if not defined $kog_id;
		
		# also, some CEGs in HMM file are not part of the 248 set 
		# and are in the bigger 458 set, so just skip these
		next unless (defined $score_cutoff{$kog});
		
		# now need to determine alignment length, first skip more lines
		<HMMSEARCH>; <HMMSEARCH>; <HMMSEARCH>; <HMMSEARCH>; <HMMSEARCH>; <HMMSEARCH>;

		# want to extract start and end coordinates of alignment
		my @start_coords;
		my @end_coords;

		# now need to see if there are multiple domains belonging to this HMM
		# just want start and end coordinates of alignments from each domain
		my $more_lines = 1;
		while($more_lines){
			$line = <HMMSEARCH>;
			if ($line =~ m/\w+/){
				my ($index, undef, $undef, undef, undef, undef, undef, undef, undef, $start, $end) = $line =~ m/\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
#				print "$index) $start-$end: $line";
				push(@start_coords, $start);
				push(@end_coords, $end);				
			} else{
				$more_lines = 0;
			}
		}
		
		# if we are here, we should have seen all possible start/end coordinates
		my $start = min(@start_coords);
		my $end   = max(@end_coords);
		my $align_length = $end - $start + 1;

#		print "start = $start, end = $end, align_len = $align_length, score = $score, score_cutoff = $score_cutoff{$kog}\n";
		if($score > $score_cutoff{$kog}){
			for my $cutoff (keys %length_cutoff){
				# for complete proteins, we want alignment length to be at least 
				# 70% of protein (from KOG) length
				my $target = ($prot_length{$kog} / 100 * $length_cutoff{$cutoff});
#				print "CUTOFF = $cutoff\tTARGET = $target\tPROT_LENGTH = $prot_length{$kog}\tLENGTH_CUTOFF $length_cutoff{$cutoff}\n";
			    if ($align_length > ($prot_length{$kog} / 100 * $length_cutoff{$cutoff})){
					$over_cutoff{$cutoff}++;
					$over_cutoff_locus{$cutoff}{$kog}++;
					$over_cutoff_group{$cutoff}{$conservation{$kog}}++;
					$over_cutoff_group_locus{$cutoff}{$conservation{$kog}}{$kog}++;
			    }			    
			}		
		} 
	}	
}

close(HMMSEARCH);


######################################################################
#
#     Computing averages and other statistics
#
######################################################################

my %average;
my %average_group;
my %locus_orth;
my %locus_orth_group;
my %percent_orth;
my %percent_orth_group;

for my $cutoff (sort keys %length_cutoff){

    if ($over_cutoff{$cutoff}){
		$average{$cutoff} = $over_cutoff{$cutoff}/ (scalar keys %{$over_cutoff_locus{$cutoff}});
    } else {
		$average{$cutoff} = 0;
    }
    
    for my $locus (keys %{$over_cutoff_locus{$cutoff}}){
		if ($over_cutoff_locus{$cutoff}{$locus} > 1) {
		    $locus_orth{$cutoff}++;
		}
    }

    if ($over_cutoff{$cutoff} && $locus_orth{$cutoff}){
		$percent_orth{$cutoff} = $locus_orth{$cutoff} / (scalar keys %{$over_cutoff_locus{$cutoff}}) * 100;
    } else { 
		$percent_orth{$cutoff} = 0;	
    }


    for my $group (sort keys %conservation_groups) {

		if ($over_cutoff_group{$cutoff}{$group}){
		    $average_group{$cutoff}{$group} =  $over_cutoff_group{$cutoff}{$group} / (scalar keys %{$over_cutoff_group_locus{$cutoff}{$group}});
		} else {
		    $average_group{$cutoff}{$group} = 0;
		}
	
		for my $locus (keys %{$over_cutoff_group_locus{$cutoff}{$group}}){
		    if ($over_cutoff_group_locus{$cutoff}{$group}{$locus} > 1) {
				$locus_orth_group{$cutoff}{$group}++;
		    }
		}

		if ($over_cutoff_group{$cutoff}{$group} && $locus_orth_group{$cutoff} && (scalar keys %{$over_cutoff_group_locus{$cutoff}{$group}})){
		    $percent_orth_group{$cutoff}{$group} = $locus_orth_group{$cutoff}{$group} / (scalar keys %{$over_cutoff_group_locus{$cutoff}{$group}}) * 100;
		} else { 
		    $percent_orth_group{$cutoff}{$group} = 0;	
		}

    }
    

}


######################################################################
#
#     Printing the output
#
######################################################################


if ($warning_version) {
    print "\n WARNING executed HMMER version has not been tested !!! \n"
}

print "\n";

print "#      Statistics of the completeness of the genome based on $total_proteins CEGs      #\n\n"; 

print  "              #Prots  %Completeness  -  #Total  Average  %Ortho \n";

print "\n";


for my $cutoff (sort keys %length_cutoff){

    if (scalar keys %{$over_cutoff_group_locus{$cutoff}} && $over_cutoff{$cutoff}){
		printf "%10s %8d %11.2f      - %5d %8.2f %9.2f\n", 
		$cutoff,
		scalar keys %{$over_cutoff_locus{$cutoff}}, 
		(scalar keys %{$over_cutoff_locus{$cutoff}})/$total_proteins * 100,
		$over_cutoff{$cutoff}, 
		$average{$cutoff},
		$percent_orth{$cutoff};
    }
    else {
		printf "%10s %8d %11.2f      - %5d %8.2f %9.2f\n",
		$cutoff, 0, 0, 0, 0, 0, 0; 
	
    }

    print "\n";

    for my $group (sort keys %conservation_groups){
		if (scalar keys %{$over_cutoff_group_locus{$cutoff}{$group}}){
			printf "   Group %s %8d %11.2f      - %5d %8.2f %9.2f\n", 
		         $group,
		         scalar keys %{$over_cutoff_group_locus{$cutoff}{$group}}, 
		         (scalar keys %{$over_cutoff_group_locus{$cutoff}{$group}}) / $conservation_groups{$group} * 100,
		         $over_cutoff_group{$cutoff}{$group}, 
			 $average_group{$cutoff}{$group},
			$percent_orth_group{$cutoff}{$group};
	    }
		else {
		    printf "   Group %s %8d %11.2f      - %5d %8.2f %9.2f\n",
		    $group, 0, 0, 0, 0, 0, 0;    
		}
    }
    print "\n";
}

print "#    These results are based on the set of genes selected by Genis Parra   #\n\n"; 
print "#    Key:                                                                  #\n"; 
print "#    Prots = number of 248 ultra-conserved CEGs present in genome          #\n"; 
print "#    %Completeness = percentage of 248 ultra-conserved CEGs present        #\n"; 
print "#    Total = total number of CEGs present including putative orthologs     #\n"; 
print "#    Average = average number of orthologs per CEG                         #\n"; 
print "#    %Ortho = percentage of detected CEGS that have more than 1 ortholog   #\n\n"; 

# print missing proteins?


my %kogs_present;

print "\n#Listing missing proteins in each category\n";
for my $cutoff (sort keys %length_cutoff){
	print "\n# Category: $cutoff \n";
	$kogs_present{$cutoff} = 248;
		
	for my $locus (sort keys %score_cutoff){
		if (!(defined $over_cutoff_locus{$cutoff}{$locus})){
			print "$locus\n";
			$kogs_present{$cutoff}--;
		}
	}	
}

# print optional warnings about how many (or few) CEGs are present?
if ($opt_v) {
	foreach my $category (sort keys %kogs_present){
		if ($category eq 'Complete' and defined $kogs_present{$category}){
			warn "NOTE: Using the more stringent cutoff criteria, $kogs_present{$category}/248 complete core eukaryotic genes (CEGs) were detected.\n";
		} elsif ($category eq 'Complete'){
			warn "NOTE: Using the more stringent cutoff criteria, 0/248 complete core eukaryotic genes (CEGs) were detected.\n";
		}
 	}
}


exit(0);

##########################################################################
##                                                                      ##
##                              THE END                                 ##
##                                                                      ##
##########################################################################
