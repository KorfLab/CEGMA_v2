###############################################################################
# geneid.pm
#
###############################################################################

use strict;
use warnings;

###############################################################################
package geneid;

# Tools to work with geneid

# Based on some tool from HMMstar.pm (Ian Korf).


sub average {
    my ($sequences) = @_;
    my $sum = 0;
    my $total = 0;
    my ($mean, $st);

    foreach my $seq (@$sequences) {
	$sum += length($seq);
	$total++;
    }
    
    $mean = $sum / $total;
    
    $sum = 0;

    foreach my $seq (@$sequences) {
	$sum += (length($seq) - $mean) * (length($seq) - $mean);

    }

    $st = sqrt($sum / $total);

    return($mean,$st); 
    
}


sub false_sites {
    my ($seq, $seed, $up, $down) = @_;
    my @false_sites;
 
    for (my $i = 0; $i < length($seq); $i++) {
	    
	my $fragment;

	if  (substr($seq, $i, length($seed)) eq $seed){
	    $fragment = substr($seq, $i-$up, length($seed) + $up + $down);
	    if (length ($fragment) == length($seed) + $up + $down) {
		push @false_sites, $fragment;
	    }
	}
    }
    return @false_sites;
};

sub filling_markov{
    my $final_length = $_[0];
    my $psuedo = $_[1];
    my %index;

    my @alpha = qw(A C G T);

    for (my $i = 0; $i < 4 ; $i++){
	print "$final_length, $alpha[$i]\n";
	filling_markov_recursively ($final_length, $psuedo, \%index, 1, $alpha[$i]);
    }

};

sub filling_markov_recursively { 
    my $final_length = $_[0];
    my $psuedo = $_[1];
    my %index = %{$_[2]};
    my $temporary_length = $_[3];
    my $string = $_[4];

    my @alpha = ('A', 'C', 'G', 'T');

    if ($final_length == $temporary_length){
	print "$string\n";
    } else {
	$temporary_length++;
	for (my $i = 0; $i < 4 ; $i++){
	    filling_markov_recursively ($final_length, $psuedo, \%index, 
					          $temporary_length, $alpha[$i]);
	}
    }
};


sub dna_word_table {
	my ($mer, $pseudo) = @_;
	
	$pseudo = 0 unless defined $pseudo;
	
	my %H;
	my $code = "";
	my @var;
	for (my $i = 0; $i < $mer; $i++) {
		my $tab = "\t" x $i;
		push @var, "\$c$i";
		$code .= "$tab foreach my \$c$i (qw(A C G T)) {\n";
	}
	$code .= "\t" x ($mer);
	my $var = join("", @var);
	$code .= "\$H{\"$var\"} = $pseudo\n";
	$code .= "}" x $mer;
		
	eval $code;
	return \%H;
}

sub log_ratio {
    my ($coding, $noncoding) = @_;
    my %compositional_model;

    for my $prefix ( sort keys %{$coding->{model}} ){
	for my $nt ( sort keys %{$coding->{model}->{$prefix}}) {
	    for my $frame ( sort keys %{$coding->{model}->{$prefix}{$nt}}){
		$compositional_model{$prefix}{$nt}{$frame} = 
		    log($coding->{model}->{$prefix}{$nt}{$frame}/$noncoding->{model}->{$prefix}{$nt}{0});

	    }
	}
    }
    return  \%compositional_model;
}

sub write_log {
    my ($compositional_model, $file) = @_; 

    open(OUT, '>', $file) or die;

    my $i = 0;

    for my $prefix ( sort keys %{$compositional_model} ){
	for my $nt ( sort keys %{$compositional_model->{$prefix}}) {
	    
	    for my $frame ( sort keys %{$compositional_model->{$prefix}{$nt}}){
		printf OUT ("%s%s %s %s %.3f\n", $prefix, $nt, $i, $frame, $compositional_model->{$prefix}{$nt}{$frame});
	    }
	    $i++;
	}
    }
    close(OUT);
		
}

sub log_ratio_pwm  {
    my ($real, $false) = @_;
    my %pwm_model;

    for my $prefix ( sort keys %{$real->{model}} ){
	for my $position ( sort keys %{$real->{model}->{$prefix}}) {
	    $pwm_model{$prefix}{$position} = 
		log($real->{model}->{$prefix}{$position}/$false->{model}->{$prefix}{$position});

	}
    }
    return  \%pwm_model;
}

sub write_pwm {
    my ($pwm_model, $file) = @_; 
    my @alpha = qw(A C G T);

    open(OUT, '>', $file) or die;
    
    my $len = scalar keys (%{$pwm_model->{A}});
    
    for (my $position = 0; $position < $len; $position++) {
	for (my $i = 0; $i < 4 ; $i++){
	    printf OUT ("%s %s %.3f\n", $position+1, $alpha[$i], $pwm_model->{$alpha[$i]}{$position});
	}
    }

    close(OUT);
		
}

###############################################################################

package geneid::SequenceModel;

sub new {
    my ($class, $name, $type, $order, $seqs, $pseudo, $frame) = @_;
    my $self = bless {
	name  => $name,
	type  => $type,
	order => $order,
    };
	
    if ($type eq 'MM') {
	$self->{model} = markov_model($order, $seqs, $pseudo, $frame);
    } elsif ($type eq 'PWM') {
	$self->{model} = pwm_model($order, $seqs, $pseudo);
    } elsif ($type eq 'FREQ') {
	$self->{model}  = frequencies_model($order, $seqs, $pseudo, $frame);    
    } else {
	die "$type unsupported at this time";
    }
	
    return $self;
}

sub write {
    my ($self, $file) = @_;
    open(OUT, '>', $file) or die;

    for my $prefix ( sort keys %{$self->{model}} ){
	for my $nt ( sort keys %{$self->{model}->{$prefix}}) {
	    for my $frame ( sort keys %{$self->{model}->{$prefix}{$nt}}){
		printf OUT ("%s %s %s %.7f\n", $prefix, $nt, $frame, $self->{model}->{$prefix}{$nt}{$frame});
	    }
	}
    }
    close(OUT);
		
}

sub frequencies_model{
    my ($order, $seqs, $pseudo, $frame) = @_;
    $pseudo = 0 unless defined $pseudo;
    my @alph = qw(A C G T);
    my $table = geneid::dna_word_table($order, "{}");
    my @total = (0, 0, 0);
    
    foreach my $prefix (keys %$table) {
	foreach my $nt (@alph) {
	    for (my $fr = 0; $fr <= $frame; $fr++){
		$table->{$prefix}{$nt}{$fr} = $pseudo;
	    }
	}
    }

	
    foreach my $seq (@$seqs) {
	my $fr=0;
	for (my $i = 0; $i < length($seq) - $order -1; $i++) {
	    my $prefix = substr($seq, $i, $order);
	    my $nt = substr($seq, $i + $order, 1);
	    next unless $prefix =~ /^[ACGT]+$/;
	    next unless $nt =~ /^[ACGT]$/;
	    $table->{$prefix}{$nt}{$fr}++;
	    $fr++; 
	     if ($fr == 3 || $frame == 0) {
		$fr = 0;
	    }
	    $total[$fr]++;
	}
    }

    foreach my $prefix (keys %$table) {
	my $zerocount = 0;

	for (my $fr = 0; $fr <= $frame; $fr++){
	    foreach my $nt (@alph) {
		if ($table->{$prefix}{$nt}{$fr} == 0) {
		    $zerocount = 1;
		}
	    }

	    if ($zerocount) {
		die "some values in Markov model with zero counts, use pseudocounts";
	    }

	    if ($total[$fr] == 0){
		die "not enough nucleotide sequence";
	    }

	    foreach my $nt (@alph) {
		$table->{$prefix}{$nt}{$fr} /= $total[$fr];
	    }
	}
    }
    return $table;
}

sub markov_model {
    my ($order, $seqs, $pseudo, $frame) = @_;
    $pseudo = 0 unless defined $pseudo;
    my @alph = qw(A C G T);
    my $table = geneid::dna_word_table($order, "{}");

    foreach my $prefix (keys %$table) {
	foreach my $nt (@alph) {
	    for (my $fr = 0; $fr <= $frame; $fr++){
		$table->{$prefix}{$nt}{$fr} = $pseudo;
	    }
	}
    }
	


    foreach my $seq (@$seqs) {
	my $fr=0;
	for (my $i = 0; $i < length($seq) - $order -1; $i++) {
	    my $prefix = substr($seq, $i, $order);
	    my $nt = substr($seq, $i + $order, 1);
	    next unless $prefix =~ /^[ACGT]+$/;
	    next unless $nt =~ /^[ACGT]$/;
	    $table->{$prefix}{$nt}{$fr}++;
	    #print "$i $prefix $nt  $fr \n";
	    $fr++; 
	    if ($fr == 3 || $frame == 0) {
		$fr = 0;
	    }
	}
    }

    foreach my $prefix (keys %$table) {
	my $zerocount = 0;
	for (my $fr = 0; $fr <= $frame; $fr++){
	    my $total = 0;
	    foreach my $nt (@alph) {
		if ($table->{$prefix}{$nt}{$fr} == 0) {
		    $zerocount = 1;
		} else {
		    $total += $table->{$prefix}{$nt}{$fr};
		}
	    }
	    if ($zerocount) {
		die "some values in Markov model with zero counts, use pseudocounts";
	    }
	
	    foreach my $nt (@alph) {
		$table->{$prefix}{$nt}{$fr} /= $total;
	    }
	}
    }
    #HMMstar::browse($table);
		
    return $table;
}

sub pwm_model{

    my ($order, $seqs, $pseudo) = @_;
    $pseudo = 0 unless defined $pseudo;
    my @alph = qw(A C G T);
    my $table = geneid::dna_word_table($order+1, "{}");

    my $len = length($seqs->[0]);
    my $total = scalar (@$seqs) + (4 * $pseudo) ;
    
       
    foreach my $prefix (keys %$table) {
	for (my $i = 0; $i < $len; $i++) {
	    $table->{$prefix}{$i} = $pseudo;
	}
    }

    foreach my $seq (@$seqs) {
	for (my $i = 0; $i < length($seq) - $order; $i++) {
	    my $prefix = substr($seq, $i, $order+1);
	    next unless $prefix =~ /^[ACGT]+$/;
	    $table->{$prefix}{$i}++;
	}
	
    }

    foreach my $prefix (keys %$table) {
	for (my $i = 0; $i < $len; $i++) {
	    if ($table->{$prefix}{$i} == 0) {
		die "some values in PWM model with zero counts, use pseudocounts";
	    }
	    
	    $table->{$prefix}{$i} /= $total;
	}
    }
		
    return $table;
}

#################################################################################



1;


=head1 Name

geneid

=head1 Description

HMMstar is an HMM-based sequence analysis software package designed for
genome analysis. This module contains a variety of functions and classes
for building traditional and generalized HMMs.

=head1 Design Principles

HMMstar is designed to be (1) useful (2) educational (3) extensible.

=head1 Examples

=cut

