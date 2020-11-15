#!/usr/bin/perl
use strict; 
use warnings;
use Getopt::Long;
use List::Util qw(max);
use Time::HiRes;

#declare vars
my $align;
my $alignfasta;
my $bps;
my $duration;
my $index;
my $loop_counter = 1;
my $maxfreq;
my $outfile = "alignment-consensus.fasta";
my $output;
my $pos;
my $row;
my $sequence_time;
my $timestart;
my $timestop;

my @alignfasta = ();
my @consensus = ();
my @seq = ();

$timestart = Time::HiRes::gettimeofday();

#Getopt
GetOptions (
    'input=s' => \$align,
);

#Verify input
if (!defined $ARGV[1] && (!defined $align || !-f $align)) {
    &print_usage("\nPlease specify an aligned file");
}

#get file contents
if (-f $align) {
    open F, '<', $align;
        @alignfasta = <F>;
    close F;
} else {&print_usage("Did you specify a file with aligned seqs?")}


#remove all elements of array that have seqis
for ($row = 0; $row <= $#alignfasta; $row++) {
    if ($alignfasta[$row] =~ m/[0-9][^GATC]*/) {
        splice(@alignfasta, $row, 1);
    }
}


#remove newlines
chomp @alignfasta;

#for each sequence in alignfasta
for (my $seqindex = 0; $seqindex <= $#alignfasta; $seqindex++) {
    @seq = ();
    @seq = split(//, $alignfasta[$seqindex]);

    $loop_counter++;

    #for each position within each sequence 
    for ($pos = 0; $pos <= $#seq; $pos++) {

        #if this position does not have a value, initialize for all bases
        $consensus[$pos]->{Afreq} = 0 unless defined $consensus[$pos]->{Afreq};
        $consensus[$pos]->{Cfreq} = 0 unless defined $consensus[$pos]->{Cfreq};
        $consensus[$pos]->{Gfreq} = 0 unless defined $consensus[$pos]->{Gfreq};
        $consensus[$pos]->{Tfreq} = 0 unless defined $consensus[$pos]->{Tfreq};
        $consensus[$pos]->{gapfreq} = 0 unless defined $consensus[$pos]->{gapfreq};

        #check contents of position
        if ($seq[$pos] eq "A") {
            $consensus[$pos]->{Afreq}++;
        }
        elsif ($seq[$pos] eq "C") {
            $consensus[$pos]->{Cfreq}++;
        }
        elsif ($seq[$pos] eq "G") {
            $consensus[$pos]->{Gfreq}++;
        }   
        elsif ($seq[$pos] eq "T") {
            $consensus[$pos]->{Tfreq}++;
        }   
        else {
            $consensus[$pos]->{gapfreq}++;
        }
    }#ends loop through $seq[$pos]
}#ends loop through $alignfasta[$seqindex]


#now to determine to position with highest frequency
for ($pos = 0; $pos <= $#consensus; $pos++) {
    $maxfreq = max ($consensus[$pos]->{Afreq}, 
                    $consensus[$pos]->{Cfreq}, 
                    $consensus[$pos]->{Gfreq}, 
                    $consensus[$pos]->{Tfreq}, 
                    $consensus[$pos]->{gapfreq});

    if ($consensus[$pos]->{Afreq} == $maxfreq) {
        $output .= "A";
    }
    elsif ($consensus[$pos]->{Cfreq} == $maxfreq) {
        $output .= "C";
    }
    elsif ($consensus[$pos]->{Gfreq} == $maxfreq) {
        $output .= "G";
    }
    elsif ($consensus[$pos]->{Tfreq} == $maxfreq) {
        $output .= "T";
    }
    elsif ($consensus[$pos]->{gapfreq} == $maxfreq) {
        $output .= "-";
    }
}


#print $output."\n";
$timestop = Time::HiRes::gettimeofday();
$duration = $timestop - $timestart;
$sequence_time = $duration / $loop_counter;
$bps = $#seq / $sequence_time;

print "This code required $sequence_time sec per sequence. \n";
print "$bps bases were calculated per second per sequence. \n";
print "total run time was $duration seconds\n";
print "$align was looked at $loop_counter times\n";

#save consensus sequence
open F, '>', $outfile or die "problem saving to file\n";
    print F ">$align\n";
    print F "$output";
close F;

###########################
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage:  --input [ALIGN.FASTA]\n";
    print "\nCheers!\n\n";
    exit;
}

