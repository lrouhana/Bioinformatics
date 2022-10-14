# Commands

## Data downloads

```
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/schizosaccharomyces_pombe/pep/Schizosaccharomyces_pombe.ASM294v2.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/protists/release-49/fasta/protists_ichthyosporea1_collection/capsaspora_owczarzaki_atcc_30864_gca_000151315/pep/Capsaspora_owczarzaki_atcc_30864_gca_000151315.C_owczarzaki_V2.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/protists/release-49/fasta/protists_choanoflagellida1_collection/salpingoeca_rosetta_gca_000188695/pep/Salpingoeca_rosetta_gca_000188695.Proterospongia_sp_ATCC50818.pep.all.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-49/fasta/amphimedon_queenslandica/pep/Amphimedon_queenslandica.Aqu1.pep.all.fa.gz

wget 'ftp://ftp.ensemblgenomes.org/pub/metazoa/release-49/fasta/nematostella_vectensis/pep/Nematostella_vectensis.ASM20922v1.pep.all.fa.gz'

wget 'ftp://ftp.ensemblgenomes.org/pub/metazoa/release-49/fasta/trichoplax_adhaerens/pep/Trichoplax_adhaerens.ASM15027v1.pep.all.fa.gz'

wget 'ftp://ftp.ensembl.org/pub/release-102/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.28.pep.all.fa.gz'

wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-49/fasta/capitella_teleta/pep/Capitella_teleta.Capitella_teleta_v1.0.pep.all.fa.gz

wget 'ftp://ftp.ensembl.org/pub/release-102/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.pep.all.fa.gz'

wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz

wget ftp://ftp.ensembl.org/pub/release-102/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz

wget https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

```

NOTE: We used the PlanMine BLAST interface (https://planmine.mpibpc.mpg.de/planmine/blast.do) to identify and download the 2 Schmidtea mediterranea CPEB proteins. These had been previously described (Rouhana et al. 2017).

<cite>Rouhana L, Tasaki J, Saberi A, Newmark PA. Genetic dissection of the planarian reproductive system through characterization of Schmidtea mediterranea CPEB homologs. Developmental biology. 2017 Jun 1;426(1):43-55.</cite>

## Identify CPEBs

```
gzip -dc *.fa.gz > all.fasta

makeblastdb -dbtype prot -in ../00-DATA/all.fasta

blastp -query Hsap_CPEB1_RBD.fa -outfmt 6 -evalue .05 -db all.fasta > blast.all.out

perl get_from_org_fa.pl > putative_cpebs_all.fa

perl isoformless.pl > isoformless_all.fa

perl rename_all.pl > rename_all.fa

perl remove_dmel_atha_drer_mmus_duplicates.pl > remove_dmel_atha_drer_mmus_duplicates.fa

mv remove_dmel_atha_drer_mmus_duplicates.fa FINAL.fa

```

## Align and run RRM tree

```

# Download HMMs from PFAM
wget -O CEBP_ZZ.hmm https://pfam.xfam.org/family/PF16366/hmm

wget -O RRM_1.hmm https://pfam.xfam.org/family/PF00076/hmm

wget -O RRM_7.hmm https://pfam.xfam.org/family/PF16367/hmm

# Build an alignment to each of the HMMs
hmm2aln.pl --hmm=CEBP_ZZ.hmm --name=CEBP_ZZ.hmm2aln --fasta=FINAL.fa --threads=20 > CEBP_ZZ.hmm2aln.fa

hmm2aln.pl --hmm=RRM_1.hmm --name=RRM1 --fasta=FINAL.fa --threads=20 > RRM_1.hmm2aln.fa 

hmm2aln.pl --hmm=RRM_7.hmm --name=RRM7 --fasta=FINAL.fa --threads=20 > RRM_7.hmm2aln.fa 

stockholm2fasta.pl rrm7.aln > rrm7.fa

stockholm2fasta.pl rrm1.aln > rrm1.fa

iqtree -nt AUTO -s CEBP_ZZ.hmm2aln.fa

perl merge_hmmsearch_outputs.pl > merge_hmmsearch_outputs.fa

mafft merge_hmmsearch_outputs.fa > rrm_merged.fa

perl -pi -e 's/ +/_/' rrm_merged.fa

iqtree -nt AUTO -s rrm_merged.fa

iqtree -nt AUTO -s CEBP_ZZ.hmm2aln.fa
```

## Identify and determine length of tails in amplicon data

```
# this script creates files with tails and files with lengths, and can create
# HTML files (for visualizing tails) if print_html() routine is uncommented
perl get_polya.pl
```


