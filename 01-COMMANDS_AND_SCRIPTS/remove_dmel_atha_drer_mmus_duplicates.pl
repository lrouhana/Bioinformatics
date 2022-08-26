#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;

our $FASTA = 'rename_all.fa';

our %GENES = ( 'AT1G17640.3' => 'AT1G17640',
               'AT1G22330.1' => 'AT1G22330',
               'AT1G22910.1' => 'AT1G22910',
               'AT1G23860.1' => 'RSZ21',
               'AT1G58470.1' => 'RBP1',
               'AT2G24590.1' => 'RSZ22a',
               'AT2G33410.1' => 'AT2G33410',
               'AT3G07810.2' => 'AT3G07810',
               'AT3G13224.2' => 'AT3G13224',
               'AT4G14300.2' => 'AT4G14300',
               'AT4G26650.1' => 'AT4G26650',
               'AT4G27000.1' => 'ATRBP45C',
               'AT4G31580.2' => 'RSZ22',
               'AT5G47620.4' => 'AT5G47620',
               'AT5G53680.1' => 'AT5G53680',
               'AT5G55550.10' => 'AT5G55550',
               'FBpp0078974' => 'Hrb27C',
               'FBpp0078975' => 'Hrb27C',
               'FBpp0078976' => 'Hrb27C',
               'FBpp0082319' => 'sqd',
               'FBpp0083776' => 'orb',
               'FBpp0083777' => 'orb',
               'FBpp0083778' => 'orb',
               'FBpp0084255' => 'msi',
               'FBpp0084256' => 'msi',
               'FBpp0099606' => 'orb',
               'FBpp0099607' => 'orb',
               'FBpp0290617' => 'Hrb27C',
               'FBpp0290618' => 'Hrb27C',
               'FBpp0290706' => 'orb',
               'FBpp0293748' => 'msi',
               'FBpp0293749' => 'msi',
               'FBpp0297873' => 'Hrb27C',
               'FBpp0297874' => 'Hrb27C',
               'FBpp0297875' => 'Hrb27C',
               'FBpp0303901' => 'orb2',
               'FBpp0303902' => 'orb2',
               'FBpp0303903' => 'orb2',
               'FBpp0303904' => 'orb2',
               'FBpp0303905' => 'orb2',
               'FBpp0305819' => 'msi',
               'FBpp0305820' => 'msi',
               'FBpp0305821' => 'msi',
               'FBpp0306565' => 'orb',
               'FBpp0423157' => 'orb2',
               'FBpp0423158' => 'orb2',
               'ENSDARP00000005968.7' => 'hnrnpdl',
               'ENSDARP00000013803.6' => 'ncl',
               'ENSDARP00000022759.7' => 'sf3b4',
               'ENSDARP00000024594.7' => 'hnrnpa1a',
               'ENSDARP00000046576.6' => 'eif4ba',
               'ENSDARP00000052510.4' => 'hnrnpa0l',
               'ENSDARP00000052511.4' => 'hnrnpa0b',
               'ENSDARP00000053266.4' => 'hnrnpa1b',
               'ENSDARP00000073682.4' => 'cpeb4a',
               'ENSDARP00000078789.5' => 'cpeb3',
               'ENSDARP00000079725.4' => 'cpeb2',
               'ENSDARP00000087189.4' => 'cpeb2',
               'ENSDARP00000089674.2' => 'cpeb1b',
               'ENSDARP00000095525.2' => 'cpeb1a',
               'ENSDARP00000108805.2' => 'hnrnpa0b',
               'ENSDARP00000108912.2' => 'hnrnpa0a',
               'ENSDARP00000109546.2' => 'ncl',
               'ENSDARP00000109818.1' => 'cpeb4a',
               'ENSDARP00000109996.1' => 'dazap1',
               'ENSDARP00000111358.1' => 'dazap1',
               'ENSDARP00000111815.2' => 'ncl',
               'ENSDARP00000112056.2' => 'eif4ba',
               'ENSDARP00000124249.1' => 'hnrnpaba',
               'ENSDARP00000124839.1' => 'hnrnpa0b',
               'ENSDARP00000125019.1' => 'hnrnpaba',
               'ENSDARP00000130727.3' => 'dazap1',
               'ENSDARP00000135390.1' => 'eif4ba',
               'ENSDARP00000137460.1' => 'cpeb4b',
               'ENSDARP00000142461.1' => 'cpeb4a',
               'ENSDARP00000144465.2' => 'cpeb4a',
               'ENSDARP00000144667.1' => 'cpeb4a',
               'ENSDARP00000144918.1' => 'hnrnpa0b',
               'ENSDARP00000146314.1' => 'dazap1',
               'ENSDARP00000146650.1' => 'cpeb4a',
               'ENSDARP00000151028.1' => 'dazap1',
               'ENSDARP00000152490.1' => 'hnrnpdl',
               'ENSDARP00000154293.1' => 'cpeb4b',
               'ENSDARP00000154503.1' => 'cpeb4a',
               'ENSDARP00000154805.1' => 'hnrnpa0l',
               'ENSDARP00000155912.1' => 'cpeb4a',
               'ENSDARP00000156574.1' => 'hnrnpdl',
               'ENSDARP00000156829.1' => 'hnrnpdl',
               'ENSMUSP00000020543.6' => 'Cpeb4',
               'ENSMUSP00000074238.3' => 'Hnrnpab',
               'ENSMUSP00000078690.4' => 'Cpeb3',
               'ENSMUSP00000084114.4' => 'Hnrnpdl',
               'ENSMUSP00000089958.5' => 'Dazap1',
               'ENSMUSP00000095936.3' => 'Cpeb1',
               'ENSMUSP00000098807.2' => 'Hnrnpab',
               'ENSMUSP00000101000.3' => 'Dazap1',
               'ENSMUSP00000101001.1' => 'Dazap1',
               'ENSMUSP00000104731.3' => 'Hnrnpab',
               'ENSMUSP00000105039.2' => 'Cpeb4',
               'ENSMUSP00000109699.3' => 'Cpeb2',
               'ENSMUSP00000109700.2' => 'Cpeb2',
               'ENSMUSP00000115656.1' => 'Cpeb4',
               'ENSMUSP00000116172.1' => 'Cpeb3',
               'ENSMUSP00000116309.1' => 'Cpeb3',
               'ENSMUSP00000116753.1' => 'Cpeb4',
               'ENSMUSP00000117497.1' => 'Dazap1',
               'ENSMUSP00000118555.1' => 'Hnrnpdl',
               'ENSMUSP00000118723.1' => 'Cpeb3',
               'ENSMUSP00000120139.1' => 'Cpeb1',
               'ENSMUSP00000120416.1' => 'Cpeb3',
               'ENSMUSP00000121005.1' => 'Hnrnpdl',
               'ENSMUSP00000121987.1' => 'Cpeb3',
               'ENSMUSP00000122442.1' => 'Cpeb3',
               'ENSMUSP00000125857.1' => 'Cpeb2',
               'ENSMUSP00000127833.1' => 'Hnrnpd',
               'ENSMUSP00000130921.2' => 'Cpeb2',
               'ENSMUSP00000136781.1' => 'Rbm31y',
               'ENSMUSP00000137079.1' => 'Cpeb1',
               'ENSMUSP00000149481.1' => 'Cpeb3',
               'ENSMUSP00000149498.1' => 'Cpeb3',
               'ENSMUSP00000150033.1' => 'Cpeb3',
               'ENSMUSP00000150396.1' => 'Cpeb3',
               'ENSMUSP00000150583.1' => 'Cpeb3',
               'ENSMUSP00000150845.1' => 'Cpeb3',
               'ENSMUSP00000151042.1' => 'Cpeb3',
               'ENSMUSP00000151095.1' => 'Cpeb3');
               
MAIN: {
    my $fp = JFR::Fasta->new($FASTA);
    my %dm_at_dr_mm = ();
    while (my $rec = $fp->get_record()) {
        $rec->{'def'} =~ m/^>([^\.]+)\.(\S+)/ or die "unexpected:$rec->{'def'}";
        my $sp = $1;
        my $gid = $2;
        if ($sp eq 'Mmus' || $sp eq 'Atha' || $sp eq 'Dmel' || $sp eq 'Drer') {
            my $gname = $GENES{$gid} or die "unexpected: $gid";
            push @{$dm_at_dr_mm{$gname}}, {'def' => $rec->{'def'},
                                           'seq' => $rec->{'seq'}};
        } else {
            print "$rec->{'def'}\n$rec->{'seq'}\n";
        }
    }

    foreach my $gname (keys %dm_at_dr_mm) {
        my $rec1 = shift @{$dm_at_dr_mm{$gname}};
        my $best = "$rec1->{'def'}\n$rec1->{'seq'}\n";
        my $len  = length($rec1->{'seq'});
        foreach my $rec (@{$dm_at_dr_mm{$gname}}) {
            if (length($rec->{'seq'}) > $len) {
                $len = length($rec->{'seq'});
                $best = "$rec->{'def'}\n$rec->{'seq'}\n";
            }
        }    
        $best =~ m/^>([^\.]+)\.(\S+)/ or die "unexpected:$best";
        $best =~ s/\..+/.$gname/;
        print $best;
    } 
}

