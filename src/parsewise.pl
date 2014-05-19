#!/usr/bin/perl

##########################################################################
#                                                                        #
#                              parsewise                                 #
#                                                                        #
#                         CEGMA pipeline code                            #
#                                                                        #
##########################################################################

use strict; use warnings; use sigtrap;

die "

# This script takes a genewise output file and produce a gff file with 
# scores based on the logical AlnBlock alignment. Genewise has to be run
# with option -gff and -alb

parsewise.pl <genewise_output>

" unless (@ARGV == 1);


my $genewise_file = $ARGV[0];

######################################################################
#
#     Security check
#
######################################################################

if (!(-e "$genewise_file" )) {
    print "File does not exist: $genewise_file \n";
    exit(1);
}

######################################################################
#
#     Reading genewise output file 
#
######################################################################


open(GENEWISE, "<", "$genewise_file") or die "Can't read from $genewise_file\n";

my $sequence_id;
my $pos1 = 1;
my $pos2;
my $score;
my @alignment_gff;
my @alignment_score;
my @alignment_frameshift;
my $MATCH_STATE = 0;
my $alignment_num = 0;
my $exon = 0;
my $total_score = 0;
my $max_alignment = 0;
my $entry = 0;
my $pos2_ex;

while (<GENEWISE>){

	# Removing some incovenient symbols
    s/[,\":[\]]/ /g;

	# Splitting the genewise output line in a vector 
	my @output_line = split;
	
	# Getting the sequence id from the gff output
	if (exists($output_line[2]) && $output_line[2] eq "match"){
		$sequence_id = $output_line[0];	
	}

	if (exists($output_line[3]) && $output_line[3] eq "LOOP_STATE"){

		if ($alignment_num > 0 || $MATCH_STATE) {
			# end of a match
			$pos2 = $output_line[5] + 1;

			if ($score < 0) {
				$score = 0;
			}
			$score = int($score*100)/100;
			$total_score += $score;
			$alignment_gff[$alignment_num][$exon] = 
			      "$sequence_id\tgenewise\tcds\t$pos1\t$pos2\t$score\t+\t.\t$sequence_id\n";
			$alignment_score[$alignment_num] = $total_score;
			if ($alignment_num > 0 && 
				$total_score > $alignment_score[$alignment_num-1]){
				$max_alignment = $alignment_num;
			} 

			#print "$sequence_id $alignment_num $total_score $max_alignment $pos1\n"; 
			
			# reseting all the alignment variables
			$MATCH_STATE = 0;
			$exon = 0;
			$alignment_num++;
			$score = 0;
			$total_score = 0;			
		} 
		
		# begining of a match
		$pos1 = $output_line[6] + 2;
		$MATCH_STATE = 1;
		# to keep for the exeptional cases
		$pos2_ex = $output_line[5] + 1;

	}
	
	if (exists($output_line[7]) && $output_line[7] =~ /5SS/){
		# begining of a match
		if ($output_line[7] =~ /PHASE_1/){
		 $pos2 = $output_line[5] + 2;
		} elsif ($output_line[7] =~ /PHASE_2/) {
			$pos2 = $output_line[5] + 3;
		} else {
			$pos2 = $output_line[5] + 1;		
		}
	}

	if (exists($output_line[7]) && $output_line[7] =~ /3SS/){
			
		# end of a match
		if ($score < 0) {
			$score = 0;
		}
		$score=int($score*100)/100;
		$alignment_gff[$alignment_num][$exon] = 
			      "$sequence_id\tgenewise\tcds\t$pos1\t$pos2\t$score\t+\t.\t$sequence_id\n";
		$exon++;
		$total_score += $score;
		$score = 0;

		# begining of a match

		if ($output_line[7] =~ /PHASE_1/){
			$pos1 = $output_line[6];
		} elsif ($output_line[7] =~ /PHASE_2/) {
			$pos1 = $output_line[6] + 1;
		} else {
			$pos1 = $output_line[6] + 2;
		}
	}
	
	
	if (exists($output_line[7]) 
			&& ($output_line[7] eq "SEQUENCE_INSERTION" 
			    || $output_line[7] eq "SEQUENCE_DELETION") ){
		$alignment_frameshift[$alignment_num]++;
	}
	
	if (exists($output_line[7]) 
			&& ($output_line[7] eq "CODON" || $output_line[7] eq "INSERT")){
		# adding scores to the current exon
		$score += $output_line[0] + 0.5;
	}
	if (exists($output_line[1]) && $output_line[1] eq "output"){

		if ($sequence_id && $entry){
			if ($alignment_frameshift[$max_alignment]){
				print "# $sequence_id WARNING  $alignment_frameshift[$max_alignment] frameshift\n";
			}			
			#print scalar(@{$alignment_gff[$max_alignment]}) ." $alignment_score[$max_alignment]\n";

			# in very exceptional cases genewise does not include a loop state at the begining
			# then we have to do the following:
			if ($alignment_num == 0 && !(exists ($alignment_gff[$max_alignment][0]))){
				$pos1 = 1;
				$pos2 = $pos2_ex;
				# end of a match
				if ($score < 0) {
					$score = 0;
				}
				$score=int($score*100)/100;
				$alignment_gff[$alignment_num][$exon] = 
			      "$sequence_id\tgenewise\tcds\t$pos1\t$pos2\t$score\t+\t.\t$sequence_id\n";			
			}						
			print @{$alignment_gff[$max_alignment]}; 
		}

		# Reseting all the variables
		$sequence_id = "";
		$pos1 = 1;
		$pos2 = "";
		$score = 0;
        @alignment_gff = ();
        @alignment_score = ();
        @alignment_frameshift = ();
        $MATCH_STATE = 0;
        $alignment_num = 0;
        $exon = 0;
        $total_score = 0;
        $max_alignment = 0;
		
		# boolean that shows if we are not in the first entry
		$entry = 1;
	}	
}

close(GENEWISE);

# printing the last entry

if ($sequence_id){
	if ($alignment_frameshift[$max_alignment]){
		print "# $sequence_id WARNING $alignment_frameshift[$max_alignment] frameshift\n";
	}

	#print scalar(@{$alignment_gff[$max_alignment]}) ." $alignment_score[$max_alignment]\n";

	# in very exceptional cases genewise does not include a loop state at the begining
	# then we have to do the following:
	if ($alignment_num == 0 && !(exists ($alignment_gff[$max_alignment][0]))){
		$pos1 = 1;
		$pos2 = $pos2_ex;
		# end of a match
		if ($score < 0) {
			$score = 0;
		}
		$score=int($score*100)/100;
		$alignment_gff[$alignment_num][$exon] = 
		"$sequence_id\tgenewise\tcds\t$pos1\t$pos2\t$score\t+\t.\t$sequence_id\n";			
	}						
	
	
	print @{$alignment_gff[$max_alignment]}; 
}

exit(0);
