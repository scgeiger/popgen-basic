#!/usr/bin/perl
use warnings;
use strict;

# Declare variables
my $infile = $ARGV[0];
(my $filename) = ($infile =~ /([^.]+)/);
my $outfile = "2ln-" . $filename . ".fasta";
my $seqID = 1;
my $row;
my $counter;

# Verify input
if (!defined $ARGV[0]) {
    &print_usage("\nNo input file was given\n"); 
}

if (-e $outfile) {
    print "Outfile cannot be overwritten.\n";
    exit;
}

# Get file contents
if (-f $infile) {
    open F, '<', $infile;
    open O, '>', $outfile;
        $counter = 0;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            if ($row =~ /^>/) {
                if ($counter == 0) {
                    print O $row . "\n";
                    $counter++;
                }
                else {
                    print O "\n" . $row . "\n";
                }
            }
            else {
                print O $row;
            }
        }
    close O;
    close F;
} else {&print_usage("Did you specify an input file?")}

print "\nSequence(s) have been adjusted.\n\n";

# Usage
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: one-line-seq.pl [fasta file] \n\n";
    print "\n\t   Please specify a fasta file as ARGV[0].
           The file can contain multiple sequences.
           Output will be saved under 2ln-[infile].\n\n";

    print "Cheers!\n\n";
    exit;
}

