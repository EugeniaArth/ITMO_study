# HW3 – Modules, configuration and containers

Nextflow pipeline that continues HW2: QC, trimming, mapping, coverage plots, and variant calling.

## Pipeline steps

1. FastQC on raw reads
2. fastp trimming
3. FastQC on trimmed reads
4. BWA index + mapping (samtools sort/index)
5. Coverage calculation and plot
6. Variant calling with **nf-core BCFTOOLS_MPILEUP** (bcftools mpileup + call)

## Project structure

```
HW3/
├── main.nf                 # entry point
├── nextflow.config         # params + profiles
├── conf/                   # base, local, cluster, container configs
├── modules/
│   ├── local/              # processes from HW2
│   └── nf-core/bcftools/   # variant calling from nf-core
├── subworkflows/
├── envs/                   # conda yml files
├── bin/
└── docker/                 # Dockerfiles
```

## Required parameters

| Parameter   | Description                    |
|------------|--------------------------------|
| `--reads`  | Paired-end reads glob pattern  |
| `--reference` | Reference FASTA             |
| `--outdir` | Output directory               |

## Optional parameters

| Parameter   | Description                    |
|------------|--------------------------------|
| `--threads` | CPU threads (default: 4)      |
| `--conda_env` | Local conda env path (local profile) |
| `--docker_registry` | Docker Hub username (default: eugeniaarth) |

## Profiles

### local
Run on your machine using a pre-built conda environment:

```bash
conda env create -f envs/hw3_all.yml -n hw3

nextflow run main.nf -profile local \
  --reads '/path/to/*_{1,2}.fastq' \
  --reference /path/to/ref.fasta \
  --conda_env ~/miniconda3/envs/hw3 \
  --outdir results
```

### cluster
SLURM cluster with separate conda yml per process:

```bash
nextflow run main.nf -profile cluster \
  --reads '/path/to/*_{1,2}.fastq' \
  --reference /path/to/ref.fasta \
  --outdir results
```

Edit `conf/cluster.config` to set your SLURM account/queue.

### docker / singularity

```bash
# build and push images first (see docker/build_and_push.sh)
nextflow run main.nf -profile docker \
  --reads '/path/to/*_{1,2}.fastq' \
  --reference /path/to/ref.fasta \
  --outdir results
```

## Docker images

| Image | Tools |
|-------|-------|
| `hw3-fastqc:1.0` | FastQC |
| `hw3-fastp:1.0` | fastp |
| `hw3-bwa-samtools:1.0` | BWA + samtools |
| `hw3-python-plot:1.0` | matplotlib plot script |
| `hw3-bcftools:1.0` | bcftools |

Build and push:

```bash
docker login
chmod +x docker/build_and_push.sh
./docker/build_and_push.sh eugeniaarth
```

## Outputs

- `fastqc_raw/`, `fastqc_trimmed/`
- `trimmed_reads/`
- `reference_index/`, `mapping/`
- `coverage/`, `coverage_plots/`
- `variants/` — VCF files from bcftools
