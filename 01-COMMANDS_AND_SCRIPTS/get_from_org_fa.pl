#!/usr/bin/perl

$|++;
use strict;
use warnings;
use Data::Dumper;

# note: the following do not have patterns in the defline
#       'Saccharomyces_cerevisiae.R64-1-1.pep.all.fa',

our $VERSION = '0.02';

#our $FILE = 'blast.1.out';
#our $OUTFILE = 'putative_cpebs.fa';
our $FILE = 'blast.all.out';
our $OUTFILE = 'putative_cpebs_all.fa';
our $DIR = '../00-DATA';
our %FA = ('Aqu' => 'Amphimedon_queenslandica.Aqu1.pep.all.fa',
           'AT' => 'Arabidopsis_thaliana.TAIR10.pep.all.fa',
           'Capte' => 'Capitella_teleta.Capitella_teleta_v1.0.pep.all.fa',
           'KJE' => 'Capsaspora_owczarzaki_atcc_30864_gca_000151315.C_owczarzaki_V2.pep.all.fa',
           'ENSDAR' => 'Danio_rerio.GRCz11.pep.all.fa',
           'FBpp' => 'Drosophila_melanogaster.BDGP6.28.pep.all.fa',
           'ENSMUS' => 'Mus_musculus.GRCm38.pep.all.fa',
           'ED' => 'Nematostella_vectensis.ASM20922v1.pep.all.fa',
           'EGD' => 'Salpingoeca_rosetta_gca_000188695.Proterospongia_sp_ATCC50818.pep.all.fa',
           'SP' => 'Schizosaccharomyces_pombe.ASM294v2.pep.all.fa',
           'mk4' => 'schmidtea_mediterranea.PRJNA12585.WBPS15.protein.fa',
           'Triad' => 'Trichoplax_adhaerens.ASM15027v1.pep.all.fa',
           'Smed'  => 'Planmine_CPEB_homologs.fa',
           'ML' => 'ML2.2.aa');

MAIN: {
    system "mv $OUTFILE 99-OLD/";
    my $ra_ids = get_data($FILE);
    my %seen = ();
    ID: foreach my $id (@{$ra_ids}) {
        foreach my $key (keys %FA) { 
            if ($id =~ m/^$key/) {
                $seen{$id} = 1;
                print "get_seq_from_fasta.pl $DIR/$FA{$key} $id >> $OUTFILE\n";
                system "get_seq_from_fasta.pl $DIR/$FA{$key} $id >> $OUTFILE";
                next ID;
            }
        }
#        print "$id\n";
    }
    foreach my $id (@{$ra_ids}) {
        next if ($seen{$id});
        print "get_seq_from_fasta.pl $DIR/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa $id >> $OUTFILE\n";
        system "get_seq_from_fasta.pl $DIR/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa $id >> $OUTFILE";
    }
}

sub get_data {
    my $file = shift;
    my %ids = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        my @fields = split /\t/, $line;
        $ids{$fields[1]}++;
    }
    my @nr_ids = keys %ids;
    return \@nr_ids;
}
