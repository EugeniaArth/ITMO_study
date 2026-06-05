# HW4 – Multi-sample pipeline

Extends HW3 with CSV input and channel operators from lecture.

## Channel operators used

| Step | Operator | What it does |
|------|----------|--------------|
| 1 | `splitCsv` + `map` | read all samples into one channel |
| 2 | `combine` | attach sample metadata to BAM files |
| 3 | `branch` | split channel by `country` field |
| 4 | `mix` | join country groups back to one channel |

## Samplesheet

`assets/samplesheet.csv`:

```csv
sample,country,fastq_1,fastq_2
SRR8071024,USA,/path/to/SRR8071024_1.fastq,/path/to/SRR8071024_2.fastq
SRR38117627,France,/path/to/SRR38117627_1.fastq,/path/to/SRR38117627_2.fastq
```

## Flow

```
CSV → map (one channel)
  → QC + mapping (all samples)
  → combine (meta + BAM)
  → branch by country
      USA    → variant calling
      France → variant calling
  → mix (one channel)
  → FILTER_VARIANTS (stub)
```

## Run

```bash
cd HW4

# test stub
nextflow run main.nf -profile local -stub-run

# full run
nextflow run main.nf -profile local \
  --samplesheet assets/samplesheet.csv \
  --reference /Users/eugenianikonorova/Documents/ITMO/data/reference/NZ_CP012026.1.fasta
```
