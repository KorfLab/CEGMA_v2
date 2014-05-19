#!/usr/bin/perl
use strict;
use warnings;
use Cegma;
use HMMstar;
use geneid;
use Getopt::Std;
our ($opt_h, $opt_v);
getopts('hv');

die "
usage: $0 [options] <gff> <fasta> <output dir>
options:
  -h  help
  -v  verbose
" unless @ARGV == 3;

my ($GFF, $DNA, $DIR) = @ARGV;

# file checks
die "$GFF file does not exist\n" if (not -e $GFF);
die "$DNA file does not exist\n" if (not -e $DNA);

my $genome = HMMstar::Genome::read($DNA, $GFF);

system("mkdir $DIR") unless -d $DIR;

# collect sequences
my (@start, @stop, @don, @acc, @coding, @intron, @inter);
my (@false_start, @false_don, @false_acc);

foreach my $contig ($genome->contigs) {

 	foreach my $neg ($contig->intergenic_sequences) {
		push @inter, $neg;
		push @inter, HMMstar::reverse_complement($neg);
		
		# Obtaining false sites (no more than 10000 needed)

		if (scalar(@false_start) < 10000) {
		    push @false_start, geneid::false_sites($neg, "ATG", 4, 2);
		}
		if (scalar(@false_don) < 10000) {
		    push @false_don, geneid::false_sites($neg, "GT",3 ,8);
		}
		if (scalar(@false_acc) < 10000) {
		    push @false_acc, geneid::false_sites($neg, "AG", 17 ,3);
		}
	}

	foreach my $cds ($contig->coding_sequences) {
	    push @start, $cds->start_site->sequence(4, 4);
		push @stop, $cds->stop_site->sequence(3, 5);

		push @coding, $cds->{transcript};
		foreach my $acc ($cds->acceptor_sites) {
			push @acc, $acc->sequence(18, 3);
		}
		foreach my $don ($cds->donor_sites) {
			push @don, $don->sequence(3, 9);
		}
		foreach my $intron ($cds->introns) {
			push @intron, $intron->sequence;
		}
	}
}


printf ("DATA COLLECTED: %s Coding sequences containing %s introns \n", scalar (@coding), scalar(@intron));

my ($mean, $st) = geneid::average(\@intron);

my $max_intron = $mean + ($st * 2) > 6000 ? 25000 : $mean + ($st * 2);

open(OUT, ">", "$DIR/intron.max") or die "Can't write to $DIR/intron.max\n";
printf OUT ("First+:Internal+                Internal+:Terminal+             30:%d block\n", $max_intron);
close (OUT);


# compositional markov models 
print "Generating intron model\n" if ($opt_v);

my $intron_initial = new geneid::SequenceModel('intron', 'FREQ', 4,
			 \@intron, 10, 0);

# $intron_initial->write("$DIR/intron.initial.5.freq");

my $intron_transition = new geneid::SequenceModel('intron', 'MM', 5,
			\@intron, 10, 0);

# $intron_transition->write("$DIR/intron.transition.5.freq");
     



print "Generating coding model\n" if ($opt_v);

my $coding_initial = new geneid::SequenceModel('coding', 'FREQ', 4,
			 \@coding, 0.25, 2);

# $coding_initial->write("$DIR/coding.initial.5.freq");

my $coding_transition = new geneid::SequenceModel('coding', 'MM', 5,
			 \@coding, 0.25, 2);

# $coding_transition->write("$DIR/coding.transition.5.freq");


my $initial_logs = geneid::log_ratio($coding_initial,$intron_initial);

my $transition_logs = geneid::log_ratio($coding_transition,$intron_transition);

geneid::write_log($initial_logs,"$DIR/coding.initial.5.logs");

geneid::write_log($transition_logs,"$DIR/coding.transition.5.logs");



# create position weight matrices
my %pwm = (
	   start => \@start,
	   acc   => \@acc,
	   don   => \@don,
);

my %false_pwm = (
	   false_start =>  \@false_start,
	   false_acc   =>  \@false_acc,
	   false_don   =>  \@false_don,
);

foreach  my $feature (keys %pwm) {
#    print  "$feature\n" if ($opt_v);
    my $real_model = new geneid::SequenceModel($feature, 'PWM', 0,
					   $pwm{$feature}, 0.25, 0);

 #   print "false_$feature\n" if ($opt_v);

    my $false_model = new geneid::SequenceModel("false_$feature", 'PWM', 0,
					   $false_pwm{"false_$feature"}, 1, 0);

    my $pwm_log = geneid::log_ratio_pwm($real_model,$false_model);
    geneid::write_pwm($pwm_log,"$DIR/$feature.logs");
	
}






__END__
