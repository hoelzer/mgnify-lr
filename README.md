# Workflow Metagenomics

![](https://img.shields.io/badge/nextflow-19.10.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)

![](https://github.com/nanozoo/wf_metagenomics/workflows/Syntax_check/badge.svg)


Maintainer: Christian

Email: christian@nanozoo.org

# Workflow

**`.join` is not working in the dag chart yet**

![chart](figures/chart.png)

# Installation

**One time installation only, to use private nextflow repos**

* create git access token [here](https://github.com/settings/tokens)
    * [x] **repo** <- click on this option
    * give it a name (e.g. nextflow) and **Generate token**
* do ``nano ~/.nextflow/scm`` and include

```java
providers {
    github {
        user = 'username'
        password = 'Personal API token'  } }
```

# Input examples

* **one** .fastq file per sample: `--nano 'sample1.fastq'`
* paired end illumina: `--illumina 'S_41_17_Cf*.R{1,2}.fastq.gz'`

# Execution example

````
nextflow run main.nf --nano '*/*.fastq' --illumina '*.R{1,2}.fastq.gz' --extra_ill illumina_reads.csv
````