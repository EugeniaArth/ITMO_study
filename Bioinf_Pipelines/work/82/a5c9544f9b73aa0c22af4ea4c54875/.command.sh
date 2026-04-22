#!/bin/bash -ue
fastp       -i SRR38117627_1.fastq       -I SRR38117627_2.fastq       -o SRR38117627_trimmed_1.fastq.gz       -O SRR38117627_trimmed_2.fastq.gz       --thread 4       --json SRR38117627.fastp.json       --html SRR38117627.fastp.html
