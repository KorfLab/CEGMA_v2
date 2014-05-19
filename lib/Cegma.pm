package Cegma;
use strict;

sub count_seqs_in_fasta_file{
	my ($file) = @_;

	my $count = `grep -c \">\" $file`;
	chomp($count);
	$count = 0 if (not defined $count);
	return($count);
}

# just check whether some software dependency exists
# e.g. tblastn, geneid

sub check_tool_exists{
	my ($tool) = @_;
	my $command = "which $tool > /dev/null;";
	my $status = system($command);

	# check status of system call
	if($status != 0){
		warn "ERROR: $tool command could not be found or is not executable\n";	
		return 0;
	} else{
		return 1;
	}
}



# simple subroutine to print a FASTA file
sub print_fasta {

    my ($entry, $file) = @_;
    my ($id, $seq) = split (/\n/,$entry);
    my ($i, $e) = (1, length($seq));
    print $file "$id\n";
    while ($i <= $e) {
		print $file substr($seq, $i - 1, 60) . "\n"; 
		$i += 60;
    }
}

# This subrutine masks with N from 0 to pos1 and from pos2 to the end
sub print_fasta_masked {
   
    my ($entry, $file, $pos1, $pos2) = @_;

    my ($id, $seq) = split (/\n/, $entry);
    my ($i, $e) = (1, length($seq));

    ($pos1 = 0)  if ($pos1 < 0);
    ($pos2 = $e) if ($pos2 <= $pos1);

	# printer header
    print $file "$id\n";
	
	# this doesn't work properly if there is no genewise hit,
	# $pos 2 ends up = 1000 which means first 1000 nt are not masked
#	print "i = $i, e = $e, pos1 = $pos1, pos2 = $pos2, L = ", length($seq),"\n";

    while ($i < $pos1) {
        print $file "N";
        print $file "\n" if (!($i % 60));
     	$i++;
    }

    while ($i <= $e && $i <= $pos2) {
		print $file  substr($seq, $i-1, 1); 
        print $file "\n" if (!($i % 60));
		$i++;
    }

    while ($i < $e) {
        print $file "N";
        print $file "\n" if (!($i % 60));
        $i++;
    }
}


sub remove_files{
	my @files = @_;
	foreach my $file (@files){
		if (-e $file){
			unlink("$file") or warn "WARNING: Can't remove $file\n"; 
		} 
	}
}


sub reverse_complement {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr[ACGTRYMKWSBDHVacgtrymkwsdbhv]
	          [TGCAYRKMSWVHDBtgcayrkmswvhdb];
	return $seq;
}


1;