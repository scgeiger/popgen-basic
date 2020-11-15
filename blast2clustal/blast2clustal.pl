#!/usr/bin/perl
use warnings;
use strict;

# Declare variables
my $infile = $ARGV[0];
(my $filename) = ($infile =~ /([^.]+)/);
my $outfile = "formatted-" . $filename . ".blast";
my $error_file = "error-" . $filename . ".blast";
my $row;
my @temp = ();
my %seqIDs = ();
my $keys;
my $length;
my $gapfreq;
my $counter;
my $test;
# Verify input

if (!defined $ARGV[0]) {
    die "\n\t Please specify a file that needs to be adjusted\n";
}

# Get file contents
if (-f $infile) {
    open F, '<', $infile;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            @temp = split(/\t/, $row);
            $keys  = ">" . $temp[0];
            if (exists($seqIDs{$keys})) { #if exists, make .1
                $counter = 1;
                $test = $keys . ".$counter";
                if (exists($seqIDs{$test})) { #if .1, make .2
                    $counter++;
                    $test = $keys . ".$counter";
                    if (exists($seqIDs{$test})) { #if .2, make .3
                        $counter++;
                        $test = $keys . ".$counter";
                        if (exists($seqIDs{$test})) { #if .3, make.4
                            $counter++;
                            $test = $keys . ".$counter";
                                if (exists($seqIDs{$test})) { #if .4, print error
                                    print "$keys has 5 versions. Error in blast2clustal.\n";
                                }
                            $keys = $keys . ".$counter"; #is .4
                        } 
                        else {$keys = $keys . ".$counter";} #is .3
                    }
                    else {$keys = $keys . ".$counter";}     #is .2 
                }
                else {$keys = $keys . ".$counter";}                  
            }
            $seqIDs{$keys} = $temp[1];
        }
    close F;
} else {&print_usage("Did you specify a file with aligned seqs?")}

open O, '>', $outfile;
    for $keys (keys %seqIDs) {
        print O $keys . "\n";
        print O $seqIDs{$keys} . "\n";
    }
close O;

#open E, '>', $error_file;
#    for $keys (%seqIDs) {
#        $length = length($seqIDs{$keys});
#        $gapfreq = $seqIDs{$keys} =~ tr/-//;
#        if (($gapfreq / $length) > 0.7) {
#            print E $keys;
#        }
#    }
#close E;

print "Blast formatted output is ready for clustal analysis as $outfile\n";
