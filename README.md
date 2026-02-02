# poolseq-snakemake

# Prep:

## Download containers:

```bash
module load apptainer
singularity pull --name multiqc.sif https://depot.galaxyproject.org/singularity/multiqc%3A1.33--pyhdfd78af_0
singularity pull --name grenedalf.sif https://depot.galaxyproject.org/singularity/grenedalf%3A0.6.3--hbefcdb2_0
singularity pull --name popoolation2.sif https://depot.galaxyproject.org/singularity/popoolation2%3A1.201--pl5321hdfd78af_0
```
## Run command example:
```bash
snakemake --profile snakeprofile --executor slurm -p --verbose --jobs 30 --use-envmodules --slurm-logdir logs/slurm --config genome=genome/E_paludinosus.fna genome_prefix=E_paludinosus pool_sizes=config/pool_sizes.tsv
```

## Work Directory Structure

```text

project-root/
├── sif/
│   ├── multiqc.sif
│   ├── grenedalf.sif
│   └── popoolation2.sif
└── poolseq-snakemake/
    ├── Snakefile
    ├── config/
    │   ├── pool_sizes.tsv
    │   └── config.yaml
    ├── genome/
    │   └── E_paludinosus.fna
    ├── results/
    ├── logs/
    │   └── slurm/
    └── snakeprofile/
