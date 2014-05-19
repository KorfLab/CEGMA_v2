###############################################################################
# HMMstar.pm
#
# ()_()
# (-.-)  Copyright 2007 Ian Korf. All rights reserved.
# (> <),
# 
###############################################################################

use strict;
use warnings;

###############################################################################
package HMMstar;
# The base class has no constructors, just general utility functions

#################
# DNA utilities #
#################

sub reverse_complement {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr[ACGTRYMKWSBDHVacgtrymkwsdbhv]
	          [TGCAYRKMSWVHDBtgcayrkmswvhdb];
	return $seq;
}

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


#########################
# Translation utilities #
#########################

# Translation uses the standard genetic code. To use other genetic codes,
# use the edit_translation() function as needed. I should make this easier.

my %Translation = (
	'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K', 'AAT' => 'N', 
	'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T', 
	'AGA' => 'R', 'AGC' => 'S', 'AGG' => 'R', 'AGT' => 'S', 
	'ATA' => 'I', 'ATC' => 'I', 'ATG' => 'M', 'ATT' => 'I', 
	'CAA' => 'Q', 'CAC' => 'H', 'CAG' => 'Q', 'CAT' => 'H', 
	'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P', 
	'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 
	'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 
	'GAA' => 'E', 'GAC' => 'D', 'GAG' => 'E', 'GAT' => 'D', 
	'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A', 
	'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G', 
	'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V', 
	'TAA' => '*', 'TAC' => 'Y', 'TAG' => '*', 'TAT' => 'Y', 
	'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 
	'TGA' => '*', 'TGC' => 'C', 'TGG' => 'W', 'TGT' => 'C', 
	'TTA' => 'L', 'TTC' => 'F', 'TTG' => 'L', 'TTT' => 'F'
);

sub edit_translation {
	my ($codon, $aa) = @_;
	die "not a codon" unless $codon =~ /^[ACGT][ACGT][ACGT]$/;
	die "not an amino acid" unless $aa =~ /^[ACDEFGHIKLMNPQRSTVWY]$/;
	$Translation{$codon} = $aa;
}

sub translate {
	my ($seq, $start) = @_;
	
	my $trans = "";
	for (my $i = $start; $i < length($seq); $i+=3) {
		my $codon = substr($seq, $i, 3);
		last if length($codon) < 3;
		if (not exists $Translation{$codon}) {$trans .= 'X'}
		else                                 {$trans .= $Translation{$codon}}
	}
	return $trans;
}

############################
# Error checking utilities #
############################

sub is_similar {
	my ($float1, $float2) = @_;
	if (abs($float1 - $float2) < 0.001) {return 1}
	return 0;
}

sub is_probability {
	my ($val) = @_;
	if ($val >= 0 and $val <= 1) {return 1}
	else                         {return 0}
}

sub is_integer {
	my ($val) = @_;
	if ($val == int($val)) {return 1}
	else                   {return 0}
}

sub browse {
	no warnings; # in this scope only
	my ($val, $level) = @_;
	$level = 0 unless defined $level;
	my $tab = '    ' x $level;
	my $ref = ref $val;
	if    (not defined $val) {print "undef\n"}
	elsif ($ref eq 'CODE' or $ref eq 'GLOB') {
		print $val, "\n";
	}
	elsif ($ref eq 'ARRAY') {
		print "$ref\n";
		$level++;
		for(my $i = 0; $i < @{$val}; $i++) {
			print "$tab    [$i] = ";
			browse($val->[$i], $level);
		}
	}
	elsif ($ref eq 'HASH' or $ref =~ /^[A-Z]/) {
		print "$ref\n";
		$level++;
		foreach my $key (sort {$a <=> $b or $a cmp $b} keys %{$val}) {
			print "$tab    $key => ";
			browse($val->{$key}, $level);
		}
	}
	elsif ($ref) {
		print $ref, "\n";
	}
	else {
		if (length($val) > 50) {
			print substr($val, 0, 40), "... (truncated for display)", "\n"
		}
		else {
			print $val, "\n";
		}
	}
}

###############################################################################
package HMMstar::DNA;

sub new {
	my ($class, $id, $seq) = @_;
	my $self = bless {
		definition => $id,
		sequence => $seq,
		'length' => CORE::length($seq),
	};
	$self->validate;
	return $self;
}

sub read {
	my ($FH) = @_;
	my $head = <$FH>;
	my ($def) = $head =~ /^>(.+)/;
	die "improper FASTA format" unless defined $def;
	my $seq = "";
	while (<$FH>) {
		chomp;
		if (/^[ACTGRYMKWSBDHVNactgrymkwsbdhvn]+$/) {$seq .= $_}
		else {die "HMMstar::DNA::read error $_"}
	}	
	return new HMMstar::DNA($def, $seq);
}

sub validate {
	my ($self) = @_;
	
	if ($self->{sequence} !~ /^[ACTGRYMKWSBDHVNactgrymkwsbdhvn]+$/) {
		die "sequence has non-DNA symbols";
	}
	if ($self->{definition} !~ /^\S+/) {
		die "sequence has no definition";
	}
}

sub definition {
	my ($self, $val) = @_;
	return $self->{definition} unless defined $val;
	$self->{definition} = $val;
	$self->validate;
}

sub sequence {
	my ($self, $val) = @_;
	return $self->{sequence} unless defined $val;
	$self->{sequence} = $val;
	$self->validate;
}

sub length {
	my ($self) = @_;
	return $self->{'length'}
}

sub reverse_complement {
	my ($self) = @_;
	my $rc = HMMstar::reverse_complement($self->{sequence});
	return new HMMstar::DNA($self->{definition}, $rc);
}

sub translate {
	my ($self, $offset) = @_;
	return HMMstar::translate($self->{sequence}, $offset);
}



###############################################################################
package HMMstar::Feature;

sub overlap {
	my ($f1, $f2) = @_;
	my ($x1, $x2) = ($f1->{start}, $f1->{end});
	my ($y1, $y2) = ($f2->{start}, $f2->{end});
	return 1 if
		($x1 >= $y1 and $x1 <= $y2) or
		($x2 >= $y1 and $x2 <= $y2);
	return 0;
}

my %Type = (
	'CDS'      => 1,
	'Intron'   => 1,
	'Exon'     => 1,
	'First'    => 1,
	'Internal' => 1,
	'Terminal' => 1,
	'Single'   => 1,
	'Donor'    => 1,
	'Acceptor' => 1,
	'Start'    => 1,
	'Stop'     => 1,
	'Inter'    => 1,
);

sub new {
	my ($class, $source, $type, $start, $end, $strand, $score, $group) = @_;
	
	my $self = bless {
		source => $source,
		type   => $type,
		start  => $start,
		end    => $end,
		strand => $strand,
		score  => $score,
		group  => $group,
	};
	$self->validate;
	return $self;
}

sub validate {
	my ($self) = @_;
	
	my $length = length($self->{source}{sequence});
				
	# check source
	die "source not HMMstar::DNA " unless ref $self->{source} eq 'HMMstar::DNA';
	
	# check type
	die "unknown type ($self->{type})" unless defined $Type{$self->{type}};
	
	# check start
	die "start not an integer" unless HMMstar::is_integer($self->{start});
	die "non-positive start coordinate ($self->{start})" if $self->{start} < 1;
	die "start out of range" if $self->{start} > $length;
	
	# check end
	die "end not an integer" unless HMMstar::is_integer($self->{end});
	die "non-positive end coordinate ($self->{end})" if $self->{end} < 1;
	die "end out of range" if $self->{end} > $length;
	
	# check strand
	die "unknown strand ($self->{strand})" unless $self->{strand} =~ /^[+-]$/;
	die "end < start" if $self->{end} < $self->{start};
	
	# score is not checked now

	# group is optional and need not be checked

}

sub source {
	my ($self, $val) = @_;
	return $self->{source} unless defined $val;
	$self->{source} = $val;
	$self->validate;
}

sub type {
	my ($self, $val) = @_;
	return $self->{type} unless defined $val;
	$self->{type} = $val;
	$self->validate;
}

sub start {
	my ($self, $val) = @_;
	return $self->{start} unless defined $val;
	$self->{start} = $val;
	$self->validate;
}

sub end {
	my ($self, $val) = @_;
	return $self->{end} unless defined $val;
	$self->{end} = $val;
	$self->validate;
}

sub strand {
	my ($self, $val) = @_;
	return $self->{strand} unless defined $val;
	$self->{strand} = $val;
	$self->validate;
}

sub score {
	my ($self, $val) = @_;
	return $self->{score} unless defined $val;
	$self->{score} = $val;
	$self->validate;
}

sub group {
	my ($self, $val) = @_;
	return $self->{group} unless defined $val;
	$self->{group} = $val;
	$self->validate;
}

sub length {
	my ($self) = @_;
	return $self->{end} - $self->{start} + 1;
}

sub sequence {
	my ($self, $up, $down) = @_;
	$up   = 0 if not defined $up;
	$down = 0 if not defined $down;
	($up, $down) = ($down, $up) if $self->{strand} eq '-';
	my $subseq = substr(
		$self->{source}{sequence},
		$self->{start} -1 - $up,
		$self->{end} - $self->{start} +1 + $up + $down
	);
	if ($self->{strand} eq '-') {
		$subseq = HMMstar::reverse_complement($subseq);
	}
	return $subseq;
}


###############################################################################
package HMMstar::CDS;

# warning conditions
my $MIN_CDS = 100;
my $MAX_CDS = 10000;
my $MIN_INTRON = 35;
my $MAX_INTRON = 100000;
my $MIN_EXON = 16;
my $MAX_EXON = 10000;
my %SPLICE_NORMAL = ("GT..AG" => 1);
my %SPLICE_WARNING = ("GC..AG" => 1, "AT..AC" => 1);

sub best_orf {
	die "best_orf not implemented yet";
}

sub new {
	my ($class, $unsorted_exons) = @_;
	
	# exons must be sorted from this point
	my @exon = sort {$a->{start} <=> $b->{start}} @$unsorted_exons;
	
	my $self = bless {
		name => $exon[0]{group},
		strand => $exon[0]{strand},
		source => $exon[0]{source},
		exons  => \@exon,
	};
	
	# check for catastrophic errors, must abort if found
	my $mixed_strand = 0;
	my $overlap = 0;
	for (my $i = 1; $i < @exon; $i++) {
		if ($exon[0]{strand} ne $exon[$i]{strand}) {$mixed_strand = 1}
		if (HMMstar::Feature::overlap($exon[$i-1], $exon[$i])) {$overlap = 1}
	}
	if ($mixed_strand or $overlap) {
		$self->{error} = 1;
		if ($mixed_strand) {warn "mixed strand in $exon[0]{group}"}
		if ($overlap) {warn "overlapping exons in $exon[0]{group}"}
		return $self;
	}
		
	# create introns
	my @intron;
	for (my $i = 1; $i < @exon; $i++) {
		my $ib = $exon[$i-1]{end} + 1;
		my $ie = $exon[$i]{start} -1;
		push @intron, new HMMstar::Feature($self->{source}, 'Intron',
			$ib, $ie, $self->{strand}, 0, $self->{name});
	}
	$self->{introns} = \@intron;
	
	# transcribe & translate
	my $tx = "";
	foreach my $exon (@{$self->{exons}}) {$tx .= $exon->sequence}
	if ($self->{strand} eq '-') {$tx = HMMstar::reverse_complement($tx)}
	$self->{transcript} = $tx;
	$self->{protein} = HMMstar::translate($tx, 0); # will validate later
	
	
	$self->validate;
	return $self;
	
}

sub name {
	my ($self) = @_;
	return $self->{name};
}

sub strand {
	my ($self) = @_;
	return $self->{strand};
}

sub validate {
	my ($self) = @_;
	
	# check for introns that are way too short
	# do a lot of checks
}

# gene features

sub introns {
	return @{shift->{introns}};
}

sub exons {
	return @{shift->{exons}};
}

sub start_site {
	my ($self) = @_;
	
	# start site defined as the 'A' of 'ATG'
	my $pos;
	if ($self->{strand} eq '+') {
		$pos = $self->{exons}[0]{start};
	} else {
		$pos = $self->{exons}[@{$self->{exons}}-1]{end};
	}
	
	my $site = new HMMstar::Feature(
		$self->{source},
		'Start',
		$pos,
		$pos,
		$self->{strand},
		0,
		$self->{group});
	return $site;
	
}

sub stop_site {
	my ($self) = @_;
	
	# stop site defined as the first 'T' in {TAA, TGA, TAG}
	my $pos;
	if ($self->{strand} eq '+') {
		$pos = $self->{exons}[@{$self->{exons}}-1]{end} -2;
	} else {
		$pos = $self->{exons}[0]{start} +2;
	}
	
	my $site = new HMMstar::Feature(
		$self->{source},
		'Stop',
		$pos,
		$pos,
		$self->{strand},
		0,
		$self->{group});
	return $site;
}

sub acceptor_sites {
	my ($self) = @_;
	
	# acceptor defined as the 'G' in the '...AG' consensus
	my @f;
	foreach my $intron (@{$self->{introns}}) {
		my $pos;
		if ($self->{strand} eq '+') {
			$pos = $intron->{end};
		} else {
			$pos = $intron->{start};
		}
		push @f, new HMMstar::Feature(
			$self->{source},
			'Acceptor',
			$pos,
			$pos,
			$self->{strand},
			0,
			$self->{group});
	}
	
	@f = reverse @f if $self->{strand} eq '-';
	return @f;
}

sub donor_sites {
	my ($self) = @_;
	
	# donor site defined as the 'G' in the 'GT...' consensus
	my @f;
	foreach my $intron (@{$self->{introns}}) {
		my $pos;
		if ($self->{strand} eq '+') {
			$pos = $intron->{start};
		} else {
			$pos = $intron->{end};
		}
		push @f, new HMMstar::Feature(
			$self->{source},
			'Donor',
			$pos,
			$pos,
			$self->{strand},
			0,
			$self->{group});
	}
	
	@f = reverse @f if $self->{strand} eq '-';
	return @f;
}


###############################################################################
package HMMstar::Contig;

sub new {
	my ($class, $dna, $features) = @_;
	my $self = bless {
		dna => $dna,
		features => $features,
	};
	return $self;
}

sub dna {return shift->{dna}}
sub features {return @{shift->{features}}}

sub intergenic_sequences {
	my ($self) = @_;
	
	my $seq = uc $self->{dna}{sequence};
	
	# sort by group
	my %group;
	foreach my $feature (@{$self->{features}}) {
		push @{$group{$feature->{group}}}, $feature;
	}
	
	foreach my $gid (keys %group) {
		my ($min, $max) = (1e20, 0);
		foreach my $feature (@{$group{$gid}}) {
			$min = $feature->{start} if $feature->{start} < $min;
			$min = $feature->{end}   if $feature->{end}   < $min;
			$max = $feature->{start} if $feature->{start} > $max;
			$max = $feature->{end}   if $feature->{end}   > $max;
		}
		substr($seq, $min -1, $max - $min +1) = lc substr($seq, $min -1, $max - $min +1);
	}
	
	my @seq = split(/[a-z]+/, $seq);
	return @seq;
}

sub coding_sequences {
	my ($self) = @_;
	
	# sort features by group name, keep only CDS features
	my %cds;
	foreach my $feature (@{$self->{features}}) {
		next unless ($feature->{type} eq 'CDS' || $feature->{type} eq 'First' || $feature->{type} eq 'Internal' || $feature->{type} eq 'Terminal' || $feature->{type} eq 'Single');
		push @{$cds{$feature->{group}}}, $feature;
	}
		
	my @cds;
	foreach my $group (keys %cds) {
		push @cds, new HMMstar::CDS($cds{$group});
	}
	
	return @cds;
}




###############################################################################
package HMMstar::Genome;

sub read {
	my ($fasta, $gff) = @_;
	my $sequences = read_fasta($fasta);
	my $features  = read_gff($gff, $sequences);
	my $self = bless {
		sequences => $sequences,
		features  => $features,
	};
	$self->validate;
	return $self;
}

sub validate {
	my ($self) = @_;
	
	# sanity check: make sure sequence identifiers in gff are in fasta
	my @unique_keys;
	foreach my $key (keys %{$self->{features}}) {
		if (not exists $self->{sequences}{$key}) {
			push @unique_keys, $key;
		}
	}
	if (@unique_keys) {
		die "Some sequence identifiers in GFF not in FASTA: @unique_keys\n";
	}
}

sub read_gff {
	my ($file, $sources) = @_;
	
	my %seq;
	open(IN, $file) or die;
	while (<IN>) {
		next if /^#/;
		my @f = split;
		push @{$seq{$f[0]}}, new HMMstar::Feature($sources->{$f[0]},
			$f[2], $f[3], $f[4], $f[6], $f[5], $f[8]);
	}
	close IN;
	
	return \%seq;
}

sub read_fasta {
	my ($file) = @_;
		
	my %seq;
	my $id;
	open(IN, $file) or die;
	while (<IN>) {
		if (/^>(.+)/) {
			$id = $1;
		} else {
			chomp;
			$seq{$id} .= $_;
		}
	}
	
	foreach my $id (keys %seq) {
		$seq{$id} = new HMMstar::DNA($id, $seq{$id});
	}
		
	return \%seq;
}

sub contigs {
	my ($self) = @_;
	
	my @contig;
	foreach my $id (sort keys %{$self->{sequences}}) {
		push @contig, new HMMstar::Contig($self->{sequences}{$id},
			$self->{features}{$id});
	}
	return @contig;
}



###############################################################################
package HMMstar::HMM;

sub read {
	my ($FH) = @_;
	
	my $self = bless {};
	
	while (<$FH>) {
		next if /^#/;
		next unless /\S/;
		last;
	}
	
	my $head = $_;
	if ($head !~ /^(HMMstar-HMM|^denada-HMM)/) {
		die "unknown format\n";
	}
	
	# read states
	my ($type, $name, $states) = split;
	$self->{type} = $type;
	$self->{name} = $name;
	$self->{states} = [];
	for (my $i = 0; $i < $states; $i++) {
		push @{$self->{states}}, HMMstar::State::read($FH);
	}
	
	# construct transition matrix
	my %tm;
	foreach my $state (@{$self->{states}}) {
		my $s1 = $state->{name};
		my $total = 0;
		foreach my $s2 (keys %{$state->{trans}}) {
			$tm{$s1}{$s2} = $state->{trans}{$s2};
			$total += $state->{trans}{$s2};
		}
		if (not HMMstar::is_similar(1, $total)) {
			HMMstar::browse($self);
			die "transitions ($total) do not sum close enough to 1.000";
		}
	}
	$self->{tmatrix} = \%tm;
		
	return $self;
}

sub write {
	my ($self, $FH) = @_;
	die "writing not yet supported";
}


###############################################################################
package HMMstar::State;

sub new {
	my ($class, $name, $init, $term, $links, $order, $duration) = @_;
	die "constructor not finished";
}

sub read {
	my ($FH) = @_;
	
	while (<$FH>) {last if /^\S/}
	my $head = $_;
	if ($head !~ /^HMMstar-State|^denada-State/) {
		die "unsupported format";
	}
	my ($type, $name, $init, $term, $links, $order, $limit) = split;
	
	my $self = bless {
		name => $name,
		init => $init,
		term => $term,
		links => $links,
		order => $order,
		limit => $limit,
	};
	
	# read trans
	my $read = 1; # lines
	my %trans;
	while (<$FH>) {
		next unless /\S/;
		my ($state, $prob) = split;
		$trans{$state} = $prob;
		last if $read++ == $links;
	}
	$self->{trans} = \%trans;
	
	# read model
	my @model;
	$read = 4 ** ($order +1); # values
	while (<$FH>) {
		next unless /\S/;
		my @val = split;
		push @model, @val;
		last if @model == $read;
	}
	$self->{model} = HMMstar::dna_word_table($order +1, 0);
	foreach my $word (sort keys %{$self->{model}}) {
		$self->{model}{$word} = shift @model;
	}

	# read range
	my @range;
	if ($limit) {
		while (<$FH>) {
			next unless /\S/;
			my @val = split;
			push @range, @val;
			last if @range == $limit
		}
	}
	$self->{range} = @range;
		
	return $self;
	
}

sub write {
}


###############################################################################
package HMMstar::SequenceModel;

sub new {
	my ($class, $name, $type, $order, $seqs, $pseudo) = @_;
	my $self = bless {
		name  => $name,
		type  => $type,
		order => $order,
	};
	
	if ($type eq 'MM') {
		$self->{model} = markov_model($order, $seqs, $pseudo);
	} elsif ($type eq 'PWM') {
		$self->{model} = pwm_model($order, $seqs, $pseudo);
	} else {
		die "$type unsupported at this time";
	}
	
	return $self;
}

sub write {
	my ($self, $file) = @_;

	open(OUT, '>', $file) or die;

	for my $prefix ( sort keys %{$self->{model}} ){
	    for my $nt ( sort keys %{$self->{model}->{$prefix}}  ) {
		printf OUT ("%s %s %.3f\n", $prefix, $nt, $self->{model}->{$prefix}{$nt});
	    }
	}
	close(OUT);
		

}

sub markov_model {
	my ($order, $seqs, $pseudo) = @_;
	$pseudo = 0 unless defined $pseudo;
	my @alph = qw(A C G T);
	my $table = HMMstar::dna_word_table($order, "{}");
	foreach my $prefix (keys %$table) {
		foreach my $nt (@alph) {
			$table->{$prefix}{$nt} = $pseudo;
		}
	}
	
	foreach my $seq (@$seqs) {
		for (my $i = 0; $i < length($seq) - $order -1; $i++) {
			my $prefix = substr($seq, $i, $order);
			my $nt = substr($seq, $i + $order + 1, 1);
			next unless $prefix =~ /^[ACGT]+$/;
			next unless $nt =~ /^[ACGT]$/;
			$table->{$prefix}{$nt}++;
		}
	}
	
	foreach my $prefix (keys %$table) {
		my $total = 0;
		my $zerocount = 0;
		foreach my $nt (@alph) {
			if ($table->{$prefix}{$nt} == 0) {
				$zerocount = 1;
			} else {
				$total += $table->{$prefix}{$nt};
			}
		}
		if ($zerocount) {
			die "some values in Markov model with zero counts, use pseudocounts";
		}
		
		foreach my $nt (@alph) {
			$table->{$prefix}{$nt} /= $total;
		}
	}
	
	
		
	return $table;
}

###############################################################################
package HMMstar::Decoder;

my $MIN_PROB  = 0;
my $MIN_SCORE = -999;

sub new {
	my ($class, $hmm, $dna) = @_;
	
	# re-map initial and terminal probabilities and transition matrix
	my %init;
	my %term;
	my %next;
	my %prev;
	foreach my $from (keys %{$hmm->{tmatrix}}) {
		foreach my $to (keys %{$hmm->{tmatrix}{$from}}) {
			if ($from eq 'INIT') {
				$init{$to} = $hmm->{tmatrix}{'INIT'}{$to};
			} elsif ($to eq 'TERM') {
				$term{$from} = $hmm->{tmatrix}{$from}{'TERM'};
			} else {
				$next{$from}{$to} = $hmm->{tmatrix}{$from}{$to};
				$prev{$to}{$from} = $hmm->{tmatrix}{$from}{$to};
			}
		}
	}
	
	# re-map states
	my %state;
	foreach my $s (@{$hmm->{states}}) {
		$state{$s->{name}} = $s;
	}
	
	my $self = bless {
		hmm     => $hmm,
		dna     => $dna,
		init    => \%init,
		term    => \%term,
		tnext   => \%next,
		tprev   => \%prev,
		state   => \%state,
		trellis => [],
		logs    => 0,
	};
	
	return $self;
}

sub decode_prob {
	my ($self) = @_;
	
	viterbi_prob($self);
	# forward
	# backward
	# posterior
}

sub decode_log {
	# viterbi
	# forward
	# backward
	# posterior
}

sub viterbi_prob {
	my ($self) = @_; 
	
	my $t = $self->{trellis};
	
	# initialization
	foreach my $s (keys %{$self->{state}}) {
		if (defined $self->{init}{$s}) {
			$t->[0]{$s}{vprob}  = $self->{init}{$s};
			$t->[0]{$s}{vtrace} = $s;
		} else {
			$t->[0]{$s}{vprob} = $MIN_PROB;
		}
	}
	
	# induction
	for (my $i = 1; $i <= $self->{dna}->length; $i++) {
		foreach my $s (keys %{$self->{state}}) {
			my ($prob, $trace) = viterbi_prob_calc($self, $i, $s);
			$t->[$i]{$s}{vprob} = $prob;
			$t->[$i]{$s}{vtrace} = $trace;
		}
	}
	
	# termination
	my $last = @$t;
	foreach my $s (keys %{$self->{state}}) {
		if (defined $self->{term}{$s}) {
			$t->[$last]{$s}{vprob} = $t->[$last-1]{$s}{vprob} * $self->{term}{$s};
		} else {
			$t->[$last]{$s}{vprob} = 0;
		}
	}
	
	# trace back
}	

sub viterbi_prob_calc {
	my ($self, $pos, $state) = @_;
	
	die "viterbi proc $pos $state\n";
	
	# emission probability
	my $emodel_name = $self->{state}{$state}{name};
	my $emodel = $self->{hmm}{emit}{$emodel_name};
	my $econtext = $self->{hmm}{context}{$emodel_name};
	
	my $nt = substr($self->{dna}{sequence}, $pos -1, 1);
	my $eprob;
	if ($econtext == 0) { 
		$eprob = $emodel->{""}{$nt};
	} elsif ($pos <= $econtext) {
		$eprob = 0.25; # if not enough context, use ambiguous value
	} else {
		my $ctx = substr($self->{dna}{sequence}, $pos - $econtext -1, $econtext);
		$eprob = $emodel->{$ctx}{$nt};
	}
	
	#print "pos = $pos, state = $state, nt = $nt, eprob = $eprob\n";
	
	# find max
	my $max_prob = 0;
	my $max_state = "";
	
	
	foreach my $from (keys %{$self->{tprev}}) {
		my $total_prob = $eprob * $self->{tprev}{$from}{$state} * 
			$self->{trellis}[$pos -1]{$from}{vprob};
		if ($total_prob > $max_prob) {
			$max_prob = $total_prob;
			$max_state = $from;
		}
	#	print "$from tprob = $total_prob\n";
	}
	#print "max: $max_state $max_prob\n";
	
	return $max_prob, $max_state;
}

sub debug {
	my ($self) = @_;
	
	my @t = @{$self->{trellis}};
	for (my $i = 0; $i < @t; $i++) {
		print "$i\t";
		if ($i == 0) {print "init\t"}
		elsif ($i == @t -1) {print "term\t"}
		else         {print substr($self->{dna}{sequence}, $i -1, 1), "\t"}
		foreach my $s (sort keys %{$t[$i]}) {
			printf "%.2g %s\t", $t[$i]{$s}{vprob}, $t[$i]{$s}{vtrace};
		}
		print "\n";
	}
}



1;


=head1 Name

HMMstar

=head1 Description

HMMstar is an HMM-based sequence analysis software package designed for
genome analysis. This module contains a variety of functions and classes
for building traditional and generalized HMMs.

=head1 Design Principles

HMMstar is designed to be (1) useful (2) educational (3) extensible.

=head1 Examples

=cut


