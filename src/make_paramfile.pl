#!/usr/bin/perl

use strict; use warnings; use sigtrap;

die "

# This script takes a bunch of files and a template and creates a
# parameter file for geneid

make_parameterfile.pl <template> <initial_markov> <transition_markov>
      <start_matrix> <acceptor_matrix> <donor_matrix> 

" unless (@ARGV == 7);

sub prt_start;
sub prt_donor;
sub prt_acceptor;


my $template          = $ARGV[0];
my $initial_markov    = $ARGV[1];
my $transition_markov = $ARGV[2];
my $start_matrix      = $ARGV[3];
my $acceptor_matrix   = $ARGV[4];
my $donor_matrix      = $ARGV[5];
my $intron_info       = $ARGV[6];

# Security check to not overwrite the results
if (!(-e "$template" )) {
    print "File does not exist: $template \n";
    exit(1);
}

if (!(-e "$initial_markov" )) {
    print "File does not exist: $initial_markov \n";
    exit(1);
}

if (!(-e "$start_matrix" )) {
    print "File does not exist: $start_matrix\n";
    exit(1);
}


open (TEMPLATE, "<", "$template") or die "Can't read from $template\n";

while (<TEMPLATE>) {
    my $new_line = $_;
    if ($new_line =~ /Start_matrix.*/){
		prt_start($start_matrix);
    } elsif ($new_line =~ /Donor_matrix.*/){
		prt_donor($donor_matrix);
    } elsif ($new_line =~ /Acceptor_matrix.*/){
		prt_acceptor($acceptor_matrix);
    } elsif ($new_line =~ /Markov_Initial_probability_matrix/){
		print;
		system("cat $initial_markov");
    } elsif ($new_line =~ /Markov_Transition_probability_matrix/){
		print;
		system("cat $transition_markov");
    } elsif ($new_line =~ /INTRAgenic connections/){
		print;
		system("cat $intron_info");
    } else {
		print;
    }
} 


close(TEMPLATE);

exit(0);

##########################################################################
##                                                                      ##
##                              THE END                                 ##
##                                                                      ##
##########################################################################
##########################################################################


##########################################################################
##                                                                      ##
##                            SUBROUTINES                               ##
##                                                                      ##
##########################################################################


sub prt_start() {

    my ($start_matrix) = @_;

    open (START, "<", "$start_matrix") or die "Can't open $start_matrix\n";

    my $bol=1;

    while (<START>){
	if (!($_ =~ /5 A .*/) && $bol) {
	    print;
	} elsif ($_ =~ /5 A .*/){
	    $bol=0;
	    prt_start_core();
	} elsif ($_ =~ /7 T .*/) {
	    $bol=1;
	}
    }
    close(START);
}

sub prt_acceptor() {

    my ($acceptor_matrix) = @_;

    open (ACCEPTOR, "<", "$acceptor_matrix") or die "can't read from $acceptor_matrix\n";

    my $bol=1;

    while (<ACCEPTOR>){
		if (!($_ =~ /18 A .*/) && $bol) {
		    print;
		} elsif ($_ =~ /18 A .*/){
		    $bol=0;
		    prt_acceptor_core();
		} elsif ($_ =~ /19 T .*/) {
		    $bol=1;
		}
    }
    close (ACCEPTOR);
}
sub prt_donor() {

    my ($donor_matrix) = @_;

    open (DONOR, "<", "$donor_matrix") or die "can't read from $donor_matrix\n";

    my $bol=1;

    while (<DONOR>){
		if (!($_ =~ /4 A .*/) && $bol) {
		    print;
		} elsif ($_ =~ /4 A .*/){
		    $bol=0;
		    prt_donor_core();
		} elsif ($_ =~ /5 T .*/) {
		    $bol=1;
		}
    }
    close (DONOR);
}

sub prt_start_core() {
    print  <<"+++EndOfStart+++";
5 A 0
5 C -9999
5 G -9999
5 T -9999
6 A -9999
6 C -9999
6 G -9999
6 T 0
7 A -9999
7 C -9999
7 G 0
7 T -9999  
+++EndOfStart+++
}

sub prt_acceptor_core() {
    print  <<"+++EndOfAcceptor+++";
18 A 0.00
18 C -9999
18 G -9999
18 T -9999
19 A -9999
19 C -9999
19 G 0.0
19 T -9999
+++EndOfAcceptor+++
}
                                                                                           
sub prt_donor_core() {
    print  <<"+++EndOfDonor+++";
4 A -9999
4 C -9999
4 G 0.000
4 T -9999
5 A -9999
5 C -10
5 G -9999
5 T 0.000
+++EndOfDonor+++
}
  
