# MGnify-LR Workflow

![](https://img.shields.io/badge/nextflow-20.01.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/uses-conda-yellow.svg)

Maintainer: Martin H&ouml;lzer

Email: hoelzer.martin@gmail.com

# Workflow

## Nanopore-only
![chart](figures/chart_ont.png)

## Hybrid
![chart](figures/chart_hybrid.png)

## Nanopore-only with decontamination
![chart](figures/chart_ont_clean.png)

# Input examples

* **one** .fastq file per sample: `--nano 'sample1.fastq'`
* paired end illumina: `--illumina 'S_41_17_Cf*.R{1,2}.fastq.gz'`

Hybrid data needs to have matching simpleNames such as ``sample1.fastq.gz`` (ONT) and ``sample1.R{1,2}.fastq.gz`` (Illumina).

# Execution example

After pulling this git:
```
nextflow run main.nf --nano test_data/ERR3662306_1_small.fastq.gz -profile local,container,docker
```

or let nextflow pull this git:

```
nextflow run hoelzer/mgnify-lr --nano ~/.nextflow/assets/hoelzer/mgnify-lr/test_data/ERR3662306_1_small.fastq.gz --output ~/mgnify-lr_output -profile local,container,docker
```

to also decontaminate your reads and final assembly just add a ``--species`` like this:

```
nextflow run hoelzer/mgnify-lr --nano ~/.nextflow/assets/hoelzer/mgnify-lr/test_data/ERR3662306_1_small.fastq.gz --output ~/mgnify-lr_output --species eco -profile local,container,docker
```

Currently supported species are:
* hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
* mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
* eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]


# Profiles
This workflow comes with some pre-defined profiles for local execution and execution on a HPC using Conda, Docker, and Singularity. You can combine different profiles to run the pipeline on different systems (local, LSF, SLURM, Cloud) using different technologies (Conda, Docker, Singularity). 

For the most _basic_ execution on a local machine using Docker use: 
```
-profile local,container,docker
```

If you want to use Conda instead of Docker:
```
-profile local,conda
```

If you want to run on a HPC w/ LSF and Singularity use:
```
-profile lsf,container,singularity
# In that case you should also adjust --cachedir, --workdir, --databases to match your HPC environment
```

