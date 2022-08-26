#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;

#our $FASTA = 'putative_cpebs.fa';
our $FASTA = 'putative_cpebs_all.fa';

MAIN: {
    my $fp = JFR::Fasta->new($FASTA);
    my %keepers = ();
    my %count = ();
    while (my $rec = $fp->get_record()) {
        my $len = length($rec->{'seq'});
        if ($rec->{'def'} =~ m/^>(\S+)\.\d+\s/) {
            my $id = $1;
$count{$id}++;
            if ($keepers{$id} && ($len > length($keepers{$id}->{'seq'}))) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            } elsif (!$keepers{$id}) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            }
        } elsif ($rec->{'def'} =~ m/^>(\S+)\.\d+$/) {
            my $id = $1;
$count{$id}++;
            if ($keepers{$id} && ($len > length($keepers{$id}->{'seq'}))) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            } elsif (!$keepers{$id}) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            }
        } elsif ($rec->{'def'} =~ m/^>(\S+)$/) {
            my $id = $1;
$count{$id}++;
            if ($keepers{$id} && ($len > length($keepers{$id}->{'seq'}))) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            } elsif (!$keepers{$id}) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            }
        } elsif ($rec->{'def'} =~ m/^>(\S+)\s/) {
            my $id = $1;
$count{$id}++;
            if ($keepers{$id} && ($len > length($keepers{$id}->{'seq'}))) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            } elsif (!$keepers{$id}) {
                $keepers{$id} = {'def' => $rec->{'def'}, 'seq' => $rec->{'seq'}};
            }
        } else {
            die "unexpected:$rec->{'def'}";
        }  
    }
    foreach my $key (keys %keepers) {
        print "$keepers{$key}->{'def'}\n$keepers{$key}->{'seq'}\n";
    }
}


