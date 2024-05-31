#!/bin/bash

mkdir genomes

cd genomes

curl -O https://ftp.ensembl.org/pub/release-112/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.gz

gunzip Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.gz

curl -O https://ftp.ensembl.org/pub/release-110/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa.gz

gunzip Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa.gz