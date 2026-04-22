#!/bin/bash -ue
cp NZ_CP012026.1.fa reference.fa
bwa index reference.fa
samtools faidx reference.fa
