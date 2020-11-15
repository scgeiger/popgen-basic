#!/usr/bin/perl
use strict;
use warnings;

# Declare Variables
my $input_file = $ARGV[0];
my @input = ();
my $row;
my $index;
my @consensus;
my $pos;
my @seq = ();
my $most_freq;
my $output;
my $outfile = "consensus.fasta";
my @seq_data = ();

# Get file contents
if (-f $input_file) {
    open F, '<', $input_file;
        @input = <F>;
        chomp @input;
    close F;
} else {&print_usage("Please specify a fasta file with alligned sequences")}

# Remove all sequence IDs
for ($row = 0; $row <= $#input; $row++) {
    if ($input[$row] =~ s/^>//) {
        next; 
    }
    else {
        push @seq_data, $input[$row];
    }
}

@input = ();

# For each sequence in alignment
for ($index = 0; $index <= $#seq_data; $index++) {
    @seq = ();
    @seq = split(//, $seq_data[$index]);

    for ($pos = 0; $pos <= $#seq; $pos++) {
#        print "at position $pos, value is $seq[$pos]\n";
        $consensus[$pos]->{$seq[$pos]} = 0 unless defined $consensus[$pos]->{$seq[$pos]};
        $consensus[$pos]->{$seq[$pos]}++;
#        foreach my $key (keys (%{$consensus[$pos]})) {
#            print "$key: $consensus[$pos]->{$key}\t";  # Data structure appears correct;
#        }
        print "\n";
    }
}

print "Sequence compositions have been assessed\n";
my $key;

foreach $pos (0..$#consensus) {
    my $largest = 0;
    $most_freq = ();
    foreach $key (keys %{$consensus[$pos]}) {
        if ($consensus[$pos]->{$key} > $largest) {
            $largest = $consensus[$pos]->{$key};
            $most_freq = $key;
        }
    }
#    print "most_freq is: $most_freq\n";
    $output .= $most_freq;
} 
        
print "Consensus sequence has been generated\n";

open F, '>', $outfile or die "problem saving to file\n";
    print F ">$input_file\n";
    print F $output;
close F;

###########################
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: make-consensus.pl [ALIGN.FASTA]\n";
    print "\nCheers!\n\n";
    exit;
}
