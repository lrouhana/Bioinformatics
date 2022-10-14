#!/usr/bin/perl

$|++;
use strict;
use warnings;
use JFR::Translate;
use Statistics::Basic;
use Storable;
use Data::Dumper;

our $VERSION = '0.12';

our $DIR = '../01-MERGE_READS';
our $LINKER = 'GGTCACCTTGATCTGAAGC';
our $LINK_N = '3' x length($LINKER);
our $KMER_LEN = 20;
our $PERCENT_MATCH_CUTOFF = 55;
our $THREEPRIME_EXTEND = 0; # experimental (setting this to 1 might be causing an issue that leads to shorter poly a tails)

our %DATA = (
'Cyclin1' => {'files' => ['C1O.aln', 'C1E.aln'],
              'fwd' => 'CCATGTAGATTCCATATCCCAAGAGG',
              'utr' => 'TATTGTTCACCCATTGTAATATAATCCTGGATTGTACAGTATTTATAAATAAATGTGTAAATAATAACTTGGTTTATTTTCTTAACATGTTAAAGGGTATCAGCTAGAGTGAGTTTTCAAAATCCTGAGTAGCAAAATGTAATTGTTTAAGTTTAATGGTCAAGCGATAACAAATAGCACAAAGACATTACAGTGTATTAGTACTTAGTAAACAAAATAGTATTTCACTTGTTTTTGATAATCACTCTAATGAAATAGCTTTTGATCGGCAATATTATTACAAAGAAAAATATGGGACAGTGACCT',
             },
'Cylcin3' => {'files' => ['C3E.aln','C32.aln'],
              'fwd' => 'GTGGCATTTAGTACTCATTCAAGGG',
              'utr' => 'TATCTTAATCAGCTTTAAAAAATAAGGTATCAGATTACCCAATTTTAGATTACCCATATTGAGAAACATATCATTTTCATCAAGAGATAACATTTTATTTAAATTTATTTTGATAATGAACAGAGCAGCTTTGATTAACAAGCAGGAATTGTATATATTGTAAATATTCAACACTATCAACAAGTAAAGTTAGCTGCATAAGGGGAAAAAACATACAATCAAAAATAAAAATAAAAACATTGAAA',
             },
    'Mos' => {'files' => ['MO.aln','ME.aln'],
              'fwd' => 'CATTTGCGCGTCATAGTCATCG',
              'utr' => 'CAATTAATTAAAGTCTTCCAACCAACATAACAAACTAAACATTTCTTGTAACGGAAGTAAAATTCAATGAAACTGAATATAATGTAGTGTATTCCTATTTGCTGTGCTACTTTATTTTGGCGATGTTATGAAAATAAACTTCTGTATGGTGCGTAATATTATCTTCGGGGTGTGTTTATACGAAGTTTCCGTGGTCACGGGGTGGTCACGGGGTCTAGACGTATCAGTAGGC',
             }
);

MAIN : {
    foreach my $key (keys %DATA) {                      # GENE
        foreach my $file (@{$DATA{$key}->{'files'}}) {  # OVARIES, EGGS
            print "\$file = $file\n";
            my $ra_kmers = kmers($DATA{$key}->{'utr'},$KMER_LEN);
            my $rev_utr = reverse($DATA{$key}->{'utr'});
            my $ra_rev_kmers = kmers($rev_utr,$KMER_LEN);
            three_prime_extend($ra_kmers,$KMER_LEN) if ($THREEPRIME_EXTEND);
#            my ($ra_seqs,$ra_masked) = get_merged_seqs("$DIR/$file",$DATA{$key},$ra_kmers);
#            Storable::store($ra_seqs,"${file}.seqs.storable");
#            Storable::store($ra_masked,"${file}.masked.storable");
            my $ra_seqs = Storable::retrieve("${file}.seqs.storable");
            my $ra_masked = Storable::retrieve("${file}.masked.storable");
            my $ra_rev_masked = Storable::retrieve("${file}.rmasked.storable");
#my $ra_rev_masked = get_rev_masked($ra_seqs,$DATA{$key},$ra_rev_kmers);
#Storable::store($ra_rev_masked,"${file}.rmasked.storable");
            my $ra_both_masked = get_both_masked($ra_masked,$ra_rev_masked);
            my $rh_data = get_data($ra_seqs,$ra_both_masked);
            my ($ra_lens,$median_tail) = get_median_tail_len($rh_data);
            print "\$median_tail = $median_tail\n";           
            print_tail_lens($file,$ra_lens);
            print_tails($file,$rh_data);
#            print_html($file,$rh_data);
        }
    }
}

sub print_tails {
    my $file = shift;
    my $rh_d = shift;
    my @lens = ();
    open OUT, ">$file.tails.fa" or die "cannot open >$file.tails.fa:$!";
    my $max_len = 0;
    my @polyas = ();
    foreach my $sq (sort {$rh_d->{$b}->[0] <=> $rh_d->{$a}->[0]} keys %{$rh_d}) {
        next unless ($rh_d->{$sq}->[1] =~ m/N{$KMER_LEN}/);
        my @mas = split /|/, $rh_d->{$sq}->[1];
        my $percent_match = get_percent_match_non_a(\@mas);
        next if ($percent_match < $PERCENT_MATCH_CUTOFF);
        my $polya = '';
        $rh_d->{$sq}->[1] =~ m/NN([^N3]+)33/;
        $polya = $1 if ($1);
        $max_len = length($polya) unless (length($polya) < $max_len);
        push @polyas, $polya;
    }
    my $num = 0;
    foreach my $polya (sort {length($b) <=> length($a)} @polyas) {
        $num++;
        next unless (length($polya) > 5);
        my $padded_polya = pad_w_dashes($polya,$max_len);
        print OUT ">$num\n$padded_polya\n";
    }
}

sub pad_w_dashes {
    my $polya = shift;
    my $max_len = shift;
    my $len = length($polya);
    return $polya if ($len == $max_len);
    my $pad = '-' x ($max_len - $len);
    my $padded = $polya . $pad;
    return $padded;
}

sub print_tail_lens {
    my $file    = shift;
    my $ra_lens = shift;
    open OUT, ">$file.tail_lens" or die "cannot open >$file.tail_lens:$!";
    foreach my $len (@{$ra_lens}) {
        print OUT "$len\n";
    }
}

sub get_median_tail_len {
    my $rh_d = shift;
    my @lens = ();
    foreach my $sq (sort {$rh_d->{$b}->[0] <=> $rh_d->{$a}->[0]} keys %{$rh_d}) {
        next unless ($rh_d->{$sq}->[1] =~ m/N{$KMER_LEN}/);
        my @mas = split /|/, $rh_d->{$sq}->[1];
        my $percent_match = get_percent_match_non_a(\@mas);
        next if ($percent_match < $PERCENT_MATCH_CUTOFF);
        push @lens, $rh_d->{$sq}->[2];
    }
    my $median = Statistics::Basic::median(@lens);
    return(\@lens,$median);
}

# combine each directional masking and mask small stretches of unmaksed between masked
sub get_both_masked {
    my $ra_m = shift;
    my $ra_r = shift;
    my @both = ();
    for (my $i = 0; $i < @{$ra_m}; $i++) {
        my @ma = split /|/, $ra_m->[$i];
        my @rm = reverse(split( /|/, $ra_r->[$i]));
        my $bseq = '';
        for (my $j = 0; $j < @ma; $j++) {
            if ($ma[$j] eq '3' || $ma[$j] eq '5' || $ma[$j] eq 'N') {
                $bseq .= $ma[$j];
            } else {
                $bseq .= $rm[$j];
            }
        }
        $bseq =~ s/NNNNN[^N]NNNNN/NNNNNNNNNNN/;
        $bseq =~ s/NNNNN[^N][^N]NNNNN/NNNNNNNNNNNN/;
        $bseq =~ s/NNNNN[^N][^N][^N]NNNNN/NNNNNNNNNNNNN/;
        $bseq =~ s/NNNNN[^N][^N][^N][^N]NNNNN/NNNNNNNNNNNNNN/;
        push @both, $bseq;
    }
    return \@both;
}

sub get_rev_masked {
    my $ra_seqs = shift;
    my $rh_d    = shift;
    my $ra_rk    = shift;
    my @rev     = ();
    foreach my $seq (@{$ra_seqs}) {
        my $rseq = reverse $seq;
        foreach my $kmer (@{$ra_rk}) {
            $rseq = sub_w_ns($rseq,$kmer,'N');
        }
        push @rev, $rseq;
    }
    return \@rev;
}

sub get_data {
    my $ra_seqs   = shift;
    my $ra_masked = shift;
    my %data      = ();
    for (my $i = 0; $i < @{$ra_seqs}; $i++) {
        my $polya_len = get_polya_len($ra_masked->[$i]);
        my $seq_len = length($ra_masked->[$i]);
        $data{$ra_seqs->[$i]} = [$seq_len,$ra_masked->[$i],$polya_len];
    }
    return \%data;
}

sub get_polya_len {
    my $seq = shift;
    my $len = 0;
    if ($seq =~ m/N([ACTG]+)3+$/) {
        $len = length($1);
    }
    return $len;
}

sub print_html {
    my $file = shift;
    my $rh_d = shift;
    open OUT, ">$file.html" or die "cannot open >$file.html:$!";
    print OUT "<h2>$file</h2><font face=courier>\n";
    my $len_match = 0;
    my $len_mismatch = 0;
    foreach my $sq (sort {$rh_d->{$b}->[0] <=> $rh_d->{$a}->[0]} keys %{$rh_d}) {
        next unless ($rh_d->{$sq}->[1] =~ m/N{$KMER_LEN}/);
        my $polya_len = $rh_d->{$sq}->[2];
        my @seq = split /|/, $sq;
        my @mas = split /|/, $rh_d->{$sq}->[1];
        my $percent_match = get_percent_match_non_a(\@mas);
        print OUT "<p>\n";
        print OUT "$percent_match" if ($percent_match < $PERCENT_MATCH_CUTOFF);
#        warn "$file: seq != mas \n" unless (scalar(@seq) == scalar(@mas));
        if (scalar(@seq) == scalar(@mas)) {
            $len_match++;
        } else {
            my $ss = join '', @seq;
            my $mm = join '', @mas;
            print "$ss\n$mm\n\n";
            $len_mismatch++;
        }
        next unless (scalar(@seq) == scalar(@mas));
        for (my $j = 0; $j < @seq; $j++) {
            if ($mas[$j] eq '5') {
                print OUT "<font color=blue>$seq[$j]</font>";
            } elsif ($mas[$j] eq '3') {
                print OUT "<font color=purple>$seq[$j]</font>";
            } elsif ($mas[$j] eq 'N') {
                print OUT "<font color=red>$seq[$j]</font>";
            } else {
                print OUT "<b>$seq[$j]</b>";
            }
        } 
        print OUT "(polya=$polya_len)</p>\n";
    }
    print "len_match: $len_match : len_mismatch: $len_mismatch\n";
}

sub get_percent_match_non_a {
    my $ra_m = shift;
    my $total = 0;
    my $match = 0;
    foreach my $nt (@{$ra_m}) {
        next if ($nt eq '3' || $nt eq '5' || $nt eq 'A' || $nt eq 'a');
        $total++;
        $match++ if ($nt eq 'N');
    }
    my $perc = ($match / $total) * 100;
    return $perc;
}

sub get_merged_seqs {
    my $file = shift;
    my $rh_d = shift;
    my $ra_k = shift;
    my @seqs = ();
    my @masked = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        chomp $line;
        next unless $line =~ m/^seq:\s+(\S+)/;
        my $seq = $1;
        my $fwd = $rh_d->{'fwd'};
        next unless ($seq =~ m/$fwd/ && $seq =~ m/$LINKER$/);
        push @seqs, $seq;
        $seq =~ s/$LINKER$/$LINK_N/;
        $seq = sub_w_ns($seq,$fwd,'5');
        foreach my $kmer (@{$ra_k}) {
            $seq = sub_w_ns($seq,$kmer,'N');
        }
        $seq = sub_w_ns($seq,$LINKER);
        push @masked, $seq;
    }
    return (\@seqs,\@masked);
}

sub sub_w_ns {
    my $seq = shift;
    my $rep = shift;
    my $sub = shift;
    my $ns  = $sub x length($rep);
    $seq =~ s/$rep/$ns/;
    return $seq;
}

sub three_prime_extend {
    my $ra_k = shift;
    my $klen = shift;
    my $last = $ra_k->[-1];
    for (my $i = 1; $i < $klen; $i++) {
        my $sub_str = 'N' x $i;
        $last =~ s/^.{$i}/$sub_str/;
        push @{$ra_k}, $last;
    }
}

sub kmers {
    my $seq = shift or return;
    my $klen = shift;
    my $len = length $seq;

    my @kmers;
    for (my $i = 0; $i + $klen <= $len; $i++) {
        push @kmers, substr($seq, $i, $klen);
    }
    return \@kmers;
}
