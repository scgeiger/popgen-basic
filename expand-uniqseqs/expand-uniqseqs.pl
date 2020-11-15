#!/usr/bin/perl
use warnings; 
use strict;

# use this to take a file of extracted uniqseqs and reexpand to unique seqs
my $in_file = $ARGV[0];
my $seq_file;
my $out_file = "expanded-". $in_file;
my @key_array = ();
my $id;
my %in_hash = ();
my $seqid;
my %seqs = ();
my $row;
my $key;

if (!defined $ARGV[0]) {
    &print_usage("\nPlease specify file to extract seqs from.");
} 

if (-f $in_file) {
    open I, '<', $in_file;
        while ($row = <I>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            if ($row =~ s/^>//) {
                $seqid = $row;
                $seqid =~ s/^\>//d;
            }
            else {
                $in_hash{$seqid} = $row;
            }
        }
    close I;
}

open O, '>', $out_file;
    foreach $key (keys %in_hash) {
        @key_array = split (/ /, $key);
        foreach $id (@key_array) {
            print O ">$id\n";
            print O "$in_hash{$key}\n";
        }
     }
close O;

sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [FILE TO SPLIT]";
    print "\nUse this when you have multiple seqids assigned to one seq\n";
    print "\nCheers!\n\n";
    exit;
}
            


