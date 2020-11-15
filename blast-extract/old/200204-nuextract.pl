#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Term::ProgressBar;

# Declaring Variables
my $blast;
my $directory;
my @blast_coords = ();
my $buffer;

my $row;
my $seq;
my $seqID;
my %output = ();
my $temp;
my $query_ID;
my $start;
my $stop;
my ($a, $b, $f, $i);
my $filename;
my @dir_files = ();
my @no_file = ();
my $revcom;
my $path;
my $length;
my $hit;
my $reversing;
my $keys;
my $outfile   = "extracted-blast.fna";
my $errorfile = "extracted-blast.err";
my $infile;
my $seqlength;
my $progress_bar;
my $max;
my %index;
my $index_file;
my @temp = ();
my $key;
my $missing_ID;
my $query_ID_2;
my $contig = "N";

# Getopt
GetOptions (
    'blast=s' => \$blast,
    'dir=s' => \$directory,
    'buffer=i' => \$buffer,
    'index=s' => \$index_file,
);

if (!defined $ARGV[1] && (!defined $blast || !defined $directory)) {
    &print_usage("\nPlease specify the required files\n");
}

if (!defined $buffer) {
    print "No buffer was defined, and so coordinates will be extracted as listed\n";
    $buffer = 0;
}

if (!defined $index_file) {
    print "No index file was indicated. If you have any GCF sequences, they may not\n be processed correctly. Please use the script gcf-index.pl\n";
}

# Make sure directory ends in a /
if ($directory !~ /\/\z/) { 
    $directory .= "/";
}

if (-f $blast) {
    open F, '<', $blast;
        @blast_coords = <F>;
    close F;
} else {&print_usage("Did you specify a blast sequence?")};

# Save sequence ID as key, all contigs as array values
if (defined $index_file) {
   open F, '<', $index_file;
        while ($row = <F>) {
            @temp = ();
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            @temp = split(/\t/, $row);
            $key = shift @temp;
            $index{$key} = [@temp];
        }
    close F;
}

if (defined $directory) {
    opendir(DIR, $directory) or die;
        @dir_files = readdir DIR;
    closedir(DIR);
    
    foreach $temp (@blast_coords) {  #for each ID/coordinate set
        $query_ID = ();
        chomp $temp;
        ($query_ID, $a, $b) = split (/\t/, $temp);
        if ($a > $b) {
            $start = $b;
            $stop = $a;
            $revcom = "Y";
        }
        elsif ($b > $a) {
            $start = $a;
            $stop = $b;
            $revcom = "N";
        }

        $start = $start - 1 - $buffer;  #coords are assumed to be indexed from 1
        $stop = $stop - 1 + $buffer;    #coords are assumed to be indexed from 1
        $length = $stop - $start;       #start and stop are indexed from 0.

        ($filename) = grep (/$query_ID/i, @dir_files);

        # if there's an abnormal filename
        if ($filename eq ".") {
            $contig = "Y";
#            print "1. $query_ID\n";
            $query_ID =~ s/ref|//g;
#            print "2. $query_ID\n";
            $query_ID =~ s/[|]//g;
            print "(126) query_ID_2 $query_ID\n";
            $query_ID_2 = $query_ID;
            $query_ID = ">" . $query_ID;
#            print "4. $query_ID\n";
            $missing_ID = ();
 
            foreach $keys (keys %index) {
#                print "looking at key $keys\n";
                foreach $i (@{$index{$keys}}) {
                    if ($i eq $query_ID) {
                        $missing_ID = $keys;
                        print "(137) found missing $missing_ID\n";
                    }
                }
            }

            if (!defined $missing_ID) {
                push (@no_file, $query_ID);
                # should I put a next here?
                next;
                # or an else for the next clause
            }

            foreach $f (@dir_files) {
                if (index ($f, $missing_ID) != -1) {
                    $filename = $f;
                }
            }
            $path = $directory . $filename;
        }

        elsif (defined $filename && $filename ne ".") {
#            print "$filename is okay\n";
            $contig = "N";
            $path = $directory . $filename;
        }
        
        print "(163) path $path\n";
#        print "checking $path\n";
        open F, '<', $path;
            $seq = ();
            while ($row = <F>) {
                chomp $row;
                if ($row =~ /^#/) {
                    next;
                }
            # above is the same for both contig Y and N
                if ($contig eq "Y") {
                    if ($row =~ s/^>//) {
#                    print "checking $row and $query_ID\n";
                        if (index($row,$query_ID_2) != -1 && !defined $seqID) {
                            $seqID = $query_ID_2;
                            $seqID =~ s/^\>//d;
                            print "(179) seqID confirmed $seqID\n";
                            $seq = ();
                        }
                        next;
                    }

                    if (defined $seqID && $row =~ s/^>//) {
                        # then stop looping to collect sequences
                        $seqID = ();
                    }

                    if (defined $seqID) {
                        if (!defined $seq) {
                            $seq = $row;
                        }
                        else {
                            $seq .= $row;
                        }
                    }
                }
                elsif ($contig eq "N" && defined $filename && $filename ne ".") {
                    if ($row =~ s/^>//) {
                        $seqID = $row;
                        $seqID =~ s/^\>//d;
                    }
                    else {$seq .= $row};
                }
            }
            if (!defined $seq) {
                print "seq isn't defined for $query_ID, see $filename\n";
            }
            $seqlength = length($seq);
            if ($start + $length > $seqlength) {
                print "The length of substring for $query_ID is greater than the sequence length.\n";
                 die;
            }

            print "(216) for $query_ID, start is $start, length is $length\n";
            $hit = substr($seq, $start, $length);
            if ($revcom eq "Y") {
                $reversing = reverse $hit;
                $reversing =~ tr/ACGTacgt/TGCAtgca/;
                $hit = $reversing;
            }
            if ($contig eq "N") {
                $output{$query_ID} = $hit;
            }
            elsif ($contig eq "Y") {
                $output{$query_ID_2} = $hit;
            }
        close F;
    }

} else {&print_usage("Did you specify a path to fna files for extraction?\n")};


# Printing
open F, '>', $outfile or die "problem saving output to file\n";
    for $keys (keys %output) {
        print F ">" . $keys . "\n";
        print F $output{$keys} . "\n";
    }
close F;

if (@no_file) {
open F, '>', $errorfile or die "problem saving error to file\n";
    print F "#sequence IDs not found in directory\n";
    foreach $i (@no_file) {
        print F $i . "\n";
    }
close F;
}
print "Extracted sequences saved successfully\n";

###########################
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: blast-extract.pl -blast [blast output] -dir [path] -buffer [numeric] -index [file]\n";
    print "Make sure the blast output is outfmt 6, sseqID, sstart, send\n";
    print "-dir requires a path for accessing the original fasta files\n";
    print "the buffer is optional and can be undefined. I've got that covered.\n";
    print "index refers to an index of filename and contig name. There's a script for that\n";
    print "\nCheers!\n\n";
    exit;
}















