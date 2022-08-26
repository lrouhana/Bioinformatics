#!/usr/bin/perl

$|++;
use strict;
use warnings;
use JFR::Fasta;
use Data::Dumper;

# note: the following do not have patterns in the defline
#       'Caenorhabditis_elegans.WBcel235.pep.all.fa',
#       'Saccharomyces_cerevisiae.R64-1-1.pep.all.fa',

#our $FILE = 'isoformless.fa';
our $FILE = 'isoformless_all.fa';
our %FA = ('Aqu' => 'Aque',
           'AT' => 'Atha',
           'Capte' => 'Ctel',
           'KJE' => 'Cowc',
           'ENSDAR' => 'Drer',
           'FBpp' => 'Dmel',
           'ENSMUS' => 'Mmus',
           'ED' => 'Nvec',
           'EGD' => 'Sros',
           'SP' => 'Spom',
           'Triad' => 'Tadh',
           'ML' => 'Mlei');
our $DIR = '../00-DATA';
our $OUTFILE = 'putative_cpebs.fa';

MAIN: {
    my $fp = JFR::Fasta->new($FILE);
    SEQ: while (my $rec = $fp->get_record()) {
        if ($rec->{'def'} =~ m/^>(\S+) .*Source:SGD/) {
            print ">Scer.$1\n$rec->{'seq'}\n";
            next SEQ;
        }
        if ($rec->{'def'} =~ m/^>Smed/) {
            print "$rec->{'def'}\n$rec->{'seq'}\n";
            next SEQ;
        }
        foreach my $key (keys %FA) { 
            if ($rec->{'def'} =~ m/^>($key\S+)/) {
                print ">$FA{$key}.$1\n$rec->{'seq'}\n";        
                next SEQ;
            }    
        }
        die "unexpected: $rec->{'def'}\n";
    }
}

