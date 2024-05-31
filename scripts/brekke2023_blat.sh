#!/bin/bash



faToTwoBit genomes/Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa \
  blat/ARS-UCD1.3.2bit




blat \
  -t=dna \
  -q=dna \
  blat/ARS-UCD1.3.2bit \
  blat/brekke2023.fasta \
  blat/brekke2023.psl 
