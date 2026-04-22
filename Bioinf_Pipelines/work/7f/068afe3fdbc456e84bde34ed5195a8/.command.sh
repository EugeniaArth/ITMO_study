#!/bin/bash -ue
bwa mem -t 10 reference.fa SRR38117627_trimmed_1.fastq.gz SRR38117627_trimmed_2.fastq.gz       | samtools sort -@ 10 -o SRR38117627.sorted.bam

samtools index SRR38117627.sorted.bam
