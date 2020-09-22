#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#USAGE:
# sample_file -avg 25 -n_content 0 -trim 50 -target 1000
# Same as:
# sample_file --avg_phred_score 25 --n_content 0 --trimmed_length 50 --target_read_count 1000


# This script takes in a fastq file
# along with a target number of reads, average quality
# cutoff, desired n_content and trimming length.  It will then find the first reads which match the
# required parameters from that file and return them without
# having to read the whole file.

# These are the configuration options
# These are all default values.
my $average_phred_score = 25;
my $n_content = 0;
my $trimmed_length = 50;
my $target_read_count = 10000;

# Gets command line arguments
GetOptions('avg_phred_score=s' => \$average_phred_score) or die "USAGE: sample_fastq_file --avg_phred_score 25 --n_content 0 --trimmed_length 50 --target_read_count 1000\n";
GetOptions('n_content=s' => \$n_content) or die "USAGE: sample_fastq_file --avg_phred_score 25 --n_content 0 --trimmed_length 50 --target_read_count 1000\n";
GetOptions('trimmed_length=s' => \$trimmed_length) or die "USAGE: sample_fastq_file --avg_phred_score 25 --n_content 0 --trimmed_length 50 --target_read_count 1000";
GetOptions('target_read_count=s' => \$target_read_count) or die "USAGE: sample_fastq_file --avg_phred_score 25 --n_content 0 --trimmed_length 50 --target_read_count 1000";

my $header;
my $seq;
my $mid;
my $quals;

my $debug = 0;

my $valid_seqs = 0;

while (<IN>) {

        $debug and sleep(1);

        $header = $_;
        $seq = <IN>;
        $mid = <IN>;
        $quals = <IN>;
        unless ($quals) {
                die "Ran out of data in the middle of a read\n";
        }

        $debug and warn "Read\n$header$seq$mid$quals\n";

        chomp $seq;
        chomp $quals;

        # Trim them if they're too long
        if (length($seq) > $trimmed_length) {
                $seq = substr($seq,0,50);
                $quals = substr($quals,0,50);
                $debug and warn "Trimmed read length\n$header$seq\n$mid$quals\n\n";
        }

        # Check if this might be colorspace data
        if ($seq =~ /\d/) {
                die "Found numbers in reads - this is probably colorspace\n";
        }


        # Count the Ns
        my $n_count = ($seq =~ tr/Nn//);
        $debug and warn "N count = $n_count\n";

        if ($n_count > $n_content) {
                $debug and warn "N content too high - skipping\n";
                next;
        }

        my $average_quality = get_average_quality($quals);

        $debug and warn "Average quality = $average_quality\n";
        if ($average_quality < $average_phred_score) {
                        $debug and warn "Quality too low - skipping\n";
        }


        print $header,$seq,"\n",$mid,$quals,"\n";
        ++$valid_seqs;

        if ($valid_seqs == $target_read_count) {
                        $debug and warn "Found enough sequences\n";
                        last;
        }
}

sub get_average_quality {
        my ($quals) = @_;
        my @quals = split(//,$quals);

        my $qual_sum = 0;

        foreach my $qual (@quals) {
                $qual_sum += ord($qual)-33;
        }

        $qual_sum /= @quals;

        return($qual_sum);
}
