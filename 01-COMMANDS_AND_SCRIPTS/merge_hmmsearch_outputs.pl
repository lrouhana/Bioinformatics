#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use Data::Dumper;

our $RRM1 = 'rrm1.fa';
our $RRM7 = 'rrm7.fa';
#our $RRM1 = 'rrm1_subset.fa';
#our $RRM7 = 'rrm7_subset.fa';
our $ORIG_FASTA = '../01-IDENTIFY_CPEBS/FINAL.fa';

MAIN: {
    my %ids = ();
    my $rh_coords1 = get_coords($RRM1,\%ids);
    my $rh_coords7 = get_coords($RRM7,\%ids);
    foreach my $id (keys %ids) {
        expand_coords($id,$rh_coords1,$rh_coords7);
    }
    foreach my $id (keys %ids) {
        print_seqs($rh_coords1,$id,$ORIG_FASTA);
    }
}

sub print_seqs {
    my $rh_c1 = shift;
    my $id    = shift;
    my $fasta = shift;
    foreach my $ra_1 (@{$rh_c1->{$id}}) {
        system "get_seq_from_fasta.pl $fasta $id $ra_1->[0] $ra_1->[1]";
    }
}

#'Drer.hnrnpa0b' => [ [ '10', '68' ],
#                     [ '101', '169' ] ]
sub expand_coords {
    my $id    = shift;
    my $rh_c1 = shift;
    my $rh_c7 = shift;
    foreach my $ra_1 (@{$rh_c1->{$id}}) {
        foreach my $ra_7 (@{$rh_c7->{$id}}) {
            next unless (@{$ra_7});
            if (overlaps($ra_1,$ra_7)) {
print "blarg\n";
                if ($ra_7->[0] < $ra_1->[0]) {
                    $ra_1->[0] = $ra_7->[0];
                }
                if ($ra_7->[1] > $ra_1->[1]) {
                    $ra_1->[1] = $ra_7->[1];
                }
                $ra_7 = [];
            }
        }
    }
    foreach my $ra_7 (@{$rh_c7->{$id}}) {
        next unless (@{$ra_7});
        push @{$rh_c1->{$id}}, [$ra_7->[0],$ra_7->[1]];
    }
}

# examples of overlap
#   [10,100] [11,99]
#   [10,100] [9,101]
#   [10,100] [9,99]
#   [10,100] [11,101]
# examples of non-overlap
#   [10,100] [1,9]
#   [10,100] [101,200]
#overlaps([10,100],[9,99]);

#overlaps([10,100],[9,99]);
sub overlaps {
    my $ra_1 = shift;
    my $ra_2 = shift;
    if ( ($ra_1->[0] <= $ra_2->[0] && $ra_1->[1] >= $ra_2->[0]) || 
         ($ra_1->[0] >= $ra_2->[0] && $ra_1->[0] <= $ra_2->[1]) )  {
        return 1;
#         print "$ra_1->[0] $ra_1->[1] -- $ra_2->[0] $ra_2->[1]\n";
#         print "OVERLAP1\n";
    } elsif ( ($ra_1->[0] < $ra_2->[0] && $ra_1->[1] < $ra_2->[0]) ||
              ($ra_1->[0] > $ra_2->[0] && $ra_1->[0] > $ra_2->[1])) {
        return 0;
#        print "$ra_1->[0] $ra_1->[1] -- $ra_2->[0] $ra_2->[1]\n";
#        print "not overlap1\n";
    } elsif ($ra_1->[0] >= $ra_2->[0] && $ra_1->[1] <= $ra_2->[1]) {
        warn "$ra_1->[0] $ra_1->[1] -- $ra_2->[0] $ra_2->[1]\n";
        die "unexpected";
    } else {
        warn "$ra_1->[0] $ra_1->[1] -- $ra_2->[0] $ra_2->[1]\n";
        die "unexpected";
    }
}

#>Drer.ncl/289-356 [subseq from] Drer.ncl
#>Drer.ncl/379-444 [subseq from] Drer.ncl
#>Drer.ncl/467-532 [subseq from] Drer.ncl
#>Drer.ncl/556-622 [subseq from] Drer.ncl
#>Aque.Aqu2.1.26043_001/32-98 [subseq from] Aque.Aqu2.1.26043_001
sub get_coords {
    my $file = shift;
    my $rh_ids = shift;
    my %coords = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next unless ($line =~ m/^>/);
        die "unexpected $line" unless ($line =~ m|^>([^/]+)/(\d+)-(\d+)|); 
        my $id = $1;
        my $c1 = $2;
        my $c2 = $3;
        push @{$coords{$id}}, [$c1,$c2];
        $rh_ids->{$id}++;
    }
    return \%coords;
}
