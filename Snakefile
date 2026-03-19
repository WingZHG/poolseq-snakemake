###########################
##    CONFIG / GLOBALS   ##
###########################

configfile: "snakeprofile/config.yaml"

pools_map = config["pools"]
short_names = list(pools_map.keys())

# Genome must be passed at runtime
GENOME_FASTA  = config.get("genome")
GENOME_PREFIX = config.get("genome_prefix", "genome")

POOL_SIZES = config.get("pool_sizes")

if POOL_SIZES is None:
    raise ValueError(
        "Missing pool_sizes.\n"
        "Run with: --config pool_sizes=<number|path/to/pool_sizes.tsv>"
    )


if GENOME_FASTA is None:
    raise ValueError(
        "Genome FASTA not provided.\n"
        "Run with: --config genome=path/to/genome.fna genome_prefix=<name>"
    )

GENOME_DIR = "results/genome_prep"

import os

def runtime_from_fastq_buckets(r1, r2):
    total_gb = (
        os.path.getsize(r1) +
        os.path.getsize(r2)
    ) / 1024**3

    if total_gb <= 60:
        return 14 * 60
    elif total_gb <= 80:
        return 18 * 60
    elif total_gb <= 100:
        return 24 * 60
    elif total_gb <= 120:
        return 30 * 60
    else:
        return 48 * 60

def mem_request_from_bam_buckets(bam):
    size_gb = os.path.getsize(bam) / 1024**3

    if size_gb <= 60:
        return 96000      # request 32 GB
    elif size_gb <= 80:
        return 128000      # request 64 GB
    elif size_gb <= 120:
        return 168000      # request 96 GB
    elif size_gb <= 160:
        return 250000     # request 128 GB
    else:
        return 250000     # request 160 GB


def threads_from_bam_size(bam):
    size_gb = os.path.getsize(bam) / 1024**3

    if size_gb <= 40:
        return 2
    elif size_gb <= 80:
        return 2
    elif size_gb <= 120:
        return 2
    elif size_gb <= 160:
        return 2
    else:
        return 2
    

def sort_mem_per_thread_gb(bam):
    size_gb = os.path.getsize(bam) / 1024**3

    if size_gb <= 40:
        return 4
    elif size_gb <= 80:
        return 6
    elif size_gb <= 120:
        return 2
    elif size_gb <= 160:
        return 2
    else:
        return 2



#snakemake --profile snakeprofile --config genome=data/genome/sticklebgenome.fna genome_prefix=sticklebgenome all

# Fallbacks for any rule that doesn't set resources explicitly
default_resources = {
    "mem_mb": 4000,
    "time": "01:00:00",
}

# Rule-class defaults (easy to tweak)
RULE_MEM_MB = {
    "genome_prep": 40000,
    "fast_qc": 24000,
    "multi_qc": 8000,
    "trim_reads": 32000,
    "trim_qc": 8000,
    "trim_multi_qc": 8000,
    "align": 96000,
    "sort": 112000,
    "rm_dupes": 124000,
    "align_stats": 8000,
    "mpileup": 8000,
    "id_indels": 48000,
    "sync": 16000,
    "rm_indels": 48000,
    "chrom_filter": 96000,
    "filter_sync": 16000,
    "rename_sync": 8000,
    "rename_sync_nf": 8000,
    "fst_genome": 8000,
    "fst_genome_nf": 8000,
    "fst_sliding": 8000,
    "diversity_genome": 8000,
    "diversity_sliding": 8000,
}

RULE_TIME = {
    "genome_prep": "2:00:00",
    "fast_qc": "4:00:00",
    "multi_qc": "1:00:00",
    "trim_reads": "4:00:00",
    "trim_qc": "6:00:00",
    "trim_multi_qc": "1:00:00",
    "align": "16:00:00",
    "sort": "8:00:00",
    "rm_dupes": "5:00:00",
    "align_stats": "1:00:00",
    "mpileup": "36:00:00",
    "id_indels": "48:00:00",
    "sync": "18:00:00",
    "rm_indels": "48:00:00",
    "chrom_filter": "12:00:00",
    "filter_sync": "36:00:00",
    "rename_sync": "1:00:00",
    "rename_sync_nf": "1:00:00",
    "fst_genome": "2:00:00",
    "fst_genome_nf": "2:00:00",
    "fst_sliding": "2:00:00",
    "diversity_genome": "2:00:00",
    "diversity_sliding": "2:00:00",
}


###########################
##        RULE ALL       ##
###########################

rule all:
  input:
    # genome prep (forces it to be scheduled)
    f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai",

    # QC reports
    "results/multi_qc/multiqc_report.html",
    "results/trim_multi_qc/trim_multiqc_report.html",

    # stats
    "results/fst_genome/fst.csv",
    "results/fst_genome_nf/fst.csv",
    "results/diversity_genome/diversity.csv",
    "results/fst_sliding/fst.csv",
    "results/diversity_sliding/diversity.csv",

     expand("results/stats/{pool}.flagstat.txt", pool=short_names)



def hhmmss_to_min(s):
    h, m, sec = map(int, s.split(":"))
    return h*60 + m + (1 if sec > 0 else 0)


###########################
##      GENOME PREP      ##
###########################

rule genome_prep:
  resources:
    mem_mb = RULE_MEM_MB["genome_prep"],
    runtime = hhmmss_to_min(RULE_TIME["genome_prep"])
  input:
    genome = GENOME_FASTA
  output:
    genome     = f"{GENOME_DIR}/{GENOME_PREFIX}.fna",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai",
    genome_bwt = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.bwt",
    genome_pac = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.pac",
    genome_ann = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.ann",
    genome_amb = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.amb",
    genome_sa  = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.sa"
  log:
    "results/logs/genome_prep/genome_prep.log"
  envmodules:
    "StdEnv/2023",
    "bwa/0.7.18",
    "samtools/1.22.1"
  threads: config["set-resources"]["genome_prep"]["threads"]
  shell:
    r"""
    mkdir -p {GENOME_DIR} results/logs/genome_prep
    cp {input.genome} {output.genome}
    bwa index -a bwtsw {output.genome} 2> {log}
    samtools faidx {output.genome} 2>> {log}
    """


###########################
##        FASTQC         ##
###########################

rule fast_qc:
  resources:
    mem_mb = RULE_MEM_MB["fast_qc"],
    runtime = hhmmss_to_min(RULE_TIME["fast_qc"])
  input:
    r1 = lambda wc: pools_map[wc.pool]["r1"],
    r2 = lambda wc: pools_map[wc.pool]["r2"]
  output:
    html1 = "results/fast_qc/{pool}_R1_fastqc.html",
    zip1  = "results/fast_qc/{pool}_R1_fastqc.zip",
    html2 = "results/fast_qc/{pool}_R2_fastqc.html",
    zip2  = "results/fast_qc/{pool}_R2_fastqc.zip"
  log:
    "results/logs/fast_qc/{pool}.log"
  envmodules:
    "StdEnv/2020",
    "fastqc/0.11.9"
  threads: config["set-resources"]["fast_qc"]["threads"]
  shell:
    r"""
    mkdir -p results/fast_qc results/logs/fast_qc
    fastqc -t {threads} -o results/fast_qc {input.r1} {input.r2} > {log} 2>&1
    """


rule multi_qc:
  resources:
    mem_mb = RULE_MEM_MB["multi_qc"],
    runtime = hhmmss_to_min(RULE_TIME["multi_qc"])
  input:
    expand("results/fast_qc/{pool}_R1_fastqc.zip", pool=short_names) +
    expand("results/fast_qc/{pool}_R2_fastqc.zip", pool=short_names)
  output:
    "results/multi_qc/multiqc_report.html"
  log:
    "results/logs/multi_qc/multiqc.log"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  threads: config["set-resources"]["multi_qc"]["threads"]
  shell:
    r"""
    mkdir -p results/multi_qc results/logs/multi_qc
    
    apptainer exec ../sif/multiqc.sif multiqc {input} -o results/multi_qc -n multiqc_report.html -f > {log} 2>&1
    """


###########################
##      TRIMMING         ##
###########################

rule trim_reads:
  resources: 
    mem_mb = RULE_MEM_MB["trim_reads"],
    runtime = hhmmss_to_min(RULE_TIME["trim_reads"])
  input:
    r1 = lambda wc: pools_map[wc.pool]["r1"],
    r2 = lambda wc: pools_map[wc.pool]["r2"]
  output:
    r1 = "results/trim_reads/{pool}_trimmed_R1.fastq.gz",
    r2 = "results/trim_reads/{pool}_trimmed_R2.fastq.gz",
    j  = "results/trim_reads/reports/{pool}.json",
    h  = "results/trim_reads/reports/{pool}.html"
  log:
    "results/logs/trim_reads/{pool}.log"
  envmodules:
    "StdEnv/2023",
    "fastp/0.24.0"
  threads: config["set-resources"]["trim_reads"]["threads"]
  shell:
    r"""
    mkdir -p results/trim_reads/reports results/logs/trim_reads
    fastp \
      -i {input.r1} -I {input.r2} \
      -o {output.r1} -O {output.r2} \
      --detect_adapter_for_pe \
      -q 20 -l 50 -g \
      -j {output.j} -h {output.h} \
      2> {log}
    """


## Trimmed QC
rule trim_qc:
  resources: 
    mem_mb = RULE_MEM_MB["trim_qc"],
    runtime = hhmmss_to_min(RULE_TIME["trim_qc"])
  input:
    r1 = "results/trim_reads/{pool}_trimmed_R1.fastq.gz",
    r2 = "results/trim_reads/{pool}_trimmed_R2.fastq.gz",
  output:
    html1 = "results/trim_qc/{pool}_trimmed_R1_fastqc.html",
    zip1  = "results/trim_qc/{pool}_trimmed_R1_fastqc.zip",
    html2 = "results/trim_qc/{pool}_trimmed_R2_fastqc.html",
    zip2  = "results/trim_qc/{pool}_trimmed_R2_fastqc.zip"
  log:
    "results/logs/trim_qc/{pool}_qc.log"
  envmodules:
    "StdEnv/2020",
    "fastqc/0.11.9"
  threads: config["set-resources"]["trim_qc"]["threads"]
  shell:
    r"""
    mkdir -p results/trim_qc results/logs/trim_qc

    fastqc -t {threads} -o results/trim_qc {input.r1} {input.r2} > {log} 2>&1
    """


## Trim Multi QC
rule trim_multi_qc:
  resources: 
    mem_mb = RULE_MEM_MB["trim_multi_qc"],
    runtime = hhmmss_to_min(RULE_TIME["trim_multi_qc"])
  input:
    expand("results/trim_qc/{pool}_trimmed_R1_fastqc.zip", pool=short_names) +
    expand("results/trim_qc/{pool}_trimmed_R2_fastqc.zip", pool=short_names)
  output:
    report = "results/trim_multi_qc/trim_multiqc_report.html"
  log:
    "results/logs/trim_multi_qc/multi_qc.log"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  threads: config["set-resources"]["trim_multi_qc"]["threads"]
  shell:
    r"""
    mkdir -p results/trim_multi_qc results/logs/trim_multi_qc

    apptainer exec ../sif/multiqc.sif multiqc {input} -o results/trim_multi_qc -n trim_multiqc_report.html -f > {log} 2>&1
    """


rule align:
  resources:
    mem_mb = RULE_MEM_MB["align"],
    runtime = lambda wc, input: runtime_from_fastq_buckets(
        input.r1,
        input.r2
    )


  input:
    r1 = "results/trim_reads/{pool}_trimmed_R1.fastq.gz",
    r2 = "results/trim_reads/{pool}_trimmed_R2.fastq.gz",
    genome = f"{GENOME_DIR}/{GENOME_PREFIX}.fna"

  output:
    bam = "results/align/{pool}.unsorted.bam"

  log:
    "results/logs/align/{pool}.log"

  envmodules:
    "bwa/0.7.18",
    "sambamba/1.0.1"

  threads: config["set-resources"]["align"]["threads"]

  params:
    bwa_threads = lambda wc, threads: max(threads - 2, 1),
    rg = "@RG\\tID:{pool}\\tSM:{pool}\\tPL:ILLUMINA"

  shell:
    r"""
    mkdir -p results/align results/logs/align

    
    echo "Runtime bucket selected: {resources.runtime} minutes" >> {log}
    echo "FASTQ sizes:" >> {log}
    du -h {input.r1} {input.r2} >> {log}

    bwa mem -t {params.bwa_threads} -R '{params.rg}' \
        {input.genome} {input.r1} {input.r2} 2>> {log} | \
    sambamba view -S -f bam -t 2 --with-header /dev/stdin -o {output.bam} 2>> {log}

    if [[ ! -s {output.bam} ]]; then
      echo "Error: Output BAM file is empty."
      exit 1
    fi
    """


###########################
##    SORT AND INDEX     ##
###########################

rule sort_index:
  resources:
    mem_mb = lambda wc, input: mem_request_from_bam_buckets(input.bam),
    runtime = hhmmss_to_min(RULE_TIME["sort"])
  threads: 
    lambda wc, input: threads_from_bam_size(input.bam)
  input:
    bam = "results/align/{pool}.unsorted.bam"
  output:
    bam = "results/align/{pool}.bam",
    bai = "results/align/{pool}.bam.bai"
  wildcard_constraints:
    pool = r"[^\.]+"
  log:
    "results/logs/sort/{pool}.log"
  envmodules:
    "sambamba/1.0.1",
    "samtools/1.22.1"
  #threads: config["set-resources"]["sort"]["threads"]

  shell:
    r"""
    mkdir -p results/logs/sort

    sambamba sort {input.bam} -t {threads} -m 2G -o {output.bam} 2>> {log}
    sambamba index -t {threads} {output.bam} 2>> {log}
    """

## sambamba sort {input.bam} -t {threads} -m 8G -o {output.bam} 2>> {log}
## samtools sort {input.bam} -@ {threads} -m 4G -o {output.bam} -O bam 2>> {log}

###########################
##    REMOVE DUPLICATES  ##
###########################

rule rm_dupes:
  resources:
    mem_mb = RULE_MEM_MB["rm_dupes"],
    runtime = hhmmss_to_min(RULE_TIME["rm_dupes"])
  input:
    bam = "results/align/{pool}.bam",
    bai = "results/align/{pool}.bam.bai"
  output:
    bam = "results/rm_dupes/{pool}.bam",
    bai = "results/rm_dupes/{pool}.bam.bai"
  log:
    "results/logs/rm_dupes/{pool}.log"
  envmodules:
    "sambamba/1.0.1"
  threads: config["set-resources"]["rm_dupes"]["threads"]
  params:
    sort_mem_gb = lambda wc, input: sort_mem_per_thread_gb(input.bam)

  shell:
    r"""
    mkdir -p results/rm_dupes results/logs/rm_dupes

    echo "BAM size: $(du -h {input.bam})" >> {log}
    echo "Threads: {threads}" >> {log}
    echo "Memory per thread: {params.sort_mem_gb}G" >> {log}
    echo "Requested memory: {resources.mem_mb} MB" >> {log}

    sambamba markdup {input.bam} \
        -t {threads} \
        -r \
        {output.bam} 2>> {log}

    sambamba index -t {threads} {output.bam} 2>> {log}
    """


###########################
##   ALIGNMENT STATS     ##
###########################

rule align_stats:
  resources:
    mem_mb = RULE_MEM_MB["align_stats"],
    runtime = hhmmss_to_min(RULE_TIME["align_stats"])
  input:
    bam = "results/align/{pool}.bam",
    bai = "results/align/{pool}.bam.bai"
  output:
    txt = "results/stats/{pool}.flagstat.txt"
  log:
    "results/logs/stats/{pool}.log"
  envmodules:
    "sambamba/1.0.1"
  threads: config["set-resources"]["align_stats"]["threads"]
  shell:
    r"""
    mkdir -p results/stats results/logs/stats

    sambamba flagstat \
      -t {threads} \
      {input.bam} \
      > {output.txt} 2>> {log}
    """


###########################
##        MPILEUP        ##
###########################

rule mpileup:
  resources:
    mem_mb = RULE_MEM_MB["mpileup"],
    runtime = hhmmss_to_min(RULE_TIME["mpileup"])
  input:
    bams = expand("results/rm_dupes/{pool}.bam", pool=short_names),
    genome = f"{GENOME_DIR}/{GENOME_PREFIX}.fna",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai"
  output:
    mpileup = "results/mpileup/all_bams.mpileup"
  log:
    "results/logs/mpileup/mpileup.log"
  envmodules:
    "samtools/1.22.1"
  params:
    max_depth = 1000,
    min_bq = 20
  threads: config["set-resources"]["mpileup"]["threads"]
  shell:
    r"""
    mkdir -p results/mpileup results/logs/mpileup

    samtools mpileup \
      -d {params.max_depth} \
      -Q {params.min_bq} \
      -B \
      -f {input.genome} \
      {input.bams} \
      -o {output.mpileup} \
      2>> {log}
    """


###########################
##       ID INDELS       ##
###########################

rule id_indels:
  resources:
    mem_mb = RULE_MEM_MB["id_indels"],
    runtime = hhmmss_to_min(RULE_TIME["id_indels"])
  input:
    mpileup = "results/mpileup/all_bams.mpileup"
  output:
    indels = "results/id_indels/indels.gtf"
  log:
    "results/logs/id_indels/id_indels.log"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  params:
    indel_window = 5,
    min_count = 2
  threads: config["set-resources"]["id_indels"]["threads"]
  shell:
    r"""
    mkdir -p results/id_indels results/logs/id_indels

    apptainer exec ../sif/popoolation2.sif identify-indel-regions.pl \
      --input {input.mpileup} \
      --output {output.indels} \
      --indel-window {params.indel_window} \
      --min-count {params.min_count} \
      2>> {log}
    """


###########################
##          SYNC         ##
###########################

rule sync:
  resources:
    mem_mb = RULE_MEM_MB["sync"],
    runtime = hhmmss_to_min(RULE_TIME["sync"])
  input:
    mpileup = "results/mpileup/all_bams.mpileup",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai"
  output:
    sync = "results/sync/sync.sync"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  log:
    "results/logs/sync/sync.log"
  params:
    #grenedalf_path = "/home/eckertl/software/grenedalf/bin/grenedalf",
    min_qual = 20,
    out_dir = "results/sync"
  threads: config["set-resources"]["sync"]["threads"]
  shell:
    r"""
    mkdir -p results/sync results/logs/sync

    apptainer exec ../sif/grenedalf.sif grenedalf sync \
      --pileup-path {input.mpileup} \
      --pileup-min-base-qual {params.min_qual} \
      --threads {threads} \
      --reference-genome-fai {input.genome_fai} \
      --out-dir {params.out_dir} \
      2>> {log}
    """


###########################
##      MASK INDELS      ##
###########################

rule rm_indels:
  resources:
    mem_mb = RULE_MEM_MB["rm_indels"],
    runtime = hhmmss_to_min(RULE_TIME["rm_indels"])
  input:
    sync = "results/sync/sync.sync",
    indels = "results/id_indels/indels.gtf"
  output:
    sync = "results/rm_indels/rm_indels.sync"
  log:
    "results/logs/rm_indels/rm_indels.log"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  threads: config["set-resources"]["rm_indels"]["threads"]
  shell:
    r"""
    mkdir -p results/rm_indels results/logs/rm_indels

    apptainer exec ../sif/popoolation2.sif filter-sync-by-gtf.pl \
      --gtf {input.indels} \
      --input {input.sync} \
      --output {output.sync} \
      2>> {log}
    """


###########################
##     SYNC FILTER       ##
###########################

rule filter_sync:
  resources:
    mem_mb = RULE_MEM_MB["filter_sync"],
    runtime = hhmmss_to_min(RULE_TIME["filter_sync"])
  input:
    sync = "results/rm_indels/rm_indels.sync"
  output:
    sync = "results/sync_filtered/sync.filtered.sync"
  log:
    "results/logs/sync_filter/sync_filter.log"
  params:
    min_cov   = lambda wc: config["sync_filter"]["min_cov"],
    max_factor= lambda wc: config["sync_filter"]["max_cov_factor"],
    min_freq  = lambda wc: config["sync_filter"]["min_freq"]
  shell:
    r"""
    mkdir -p results/sync_filtered results/logs/sync_filter

    echo "Calculating per-pool mean coverage..." > {log}

    ##################################################
    # compute mean coverage per pool column
    ##################################################
    awk -F'\t' '
    {{
      for(i=4;i<=NF;i++)
      {{
        split($i,a,":")
        cov=a[1]+a[2]+a[3]+a[4]
        sum[i]+=cov
        n[i]++
      }}
    }}
    END{{
      for(i=4;i<=NF;i++)
        print sum[i]/n[i]
    }}
    ' {input.sync} > results/sync_filtered/pool_means.txt

    ##################################################
    # convert means → per-pool max coverage thresholds
    ##################################################
    awk -v f={params.max_factor} '{{{{print int($1*f)}}}}' \
    results/sync_filtered/pool_means.txt \
    > results/sync_filtered/pool_max.txt

    echo "Per-pool mean coverages:" >> {log}
    cat results/sync_filtered/pool_means.txt >> {log}

    echo "Per-pool max thresholds:" >> {log}
    cat results/sync_filtered/pool_max.txt >> {log}

    ##################################################
    # filtering
    ##################################################
    awk -F'\t' \
        -v mincov={params.min_cov} \
        -v minfreq={params.min_freq} \
        -v maxfile="results/sync_filtered/pool_max.txt" '
    BEGIN{{
      idx=4
      while((getline line < maxfile)>0){{
        maxcov[idx]=line
        idx++
      }}
    }}
    {{
      validPools=0
      polymorphic=0

      for(i=4;i<=NF;i++)
      {{
        split($i,a,":")
        cov=a[1]+a[2]+a[3]+a[4]

        if(cov < mincov || cov > maxcov[i])
          continue

        validPools++

        max=a[1]
        if(a[2]>max) max=a[2]
        if(a[3]>max) max=a[3]
        if(a[4]>max) max=a[4]

        if((cov-max)/cov >= minfreq)
          polymorphic=1
      }}

      if((validPools >= (0.8*(NF-3))) && polymorphic)
        print
    }}
    ' {input.sync} > {output.sync} 2>> {log}

    echo "Filtering complete." >> {log}
    """



###########################
##       RENAME SYNC     ##
###########################

rule rename_sync:
  resources:
    mem_mb = RULE_MEM_MB["rename_sync"],
    runtime = hhmmss_to_min(RULE_TIME["rename_sync"])
  input:
    sync = "results/sync_filtered/sync.filtered.sync",
    bams = expand("results/rm_dupes/{pool}.bam", pool=short_names)
  output:
    sync = "results/rename_sync/sync_final.sync"
  log:
    "results/logs/rename_sync/rename_sync.log"
  threads: config["set-resources"]["rename_sync"]["threads"]
  shell:
      r"""
      mkdir -p results/rename_sync results/logs/rename_sync

      new_headers=$(printf "%s\n" \
        BA_PTL BA_KAH BA_BWE BA_BUN BA_RSW BA_MIK-UP BA_SEB BA_NJR BA_ISW BA_BIH \
        BA_ISE BN_NJR BN_MIK-UP BN_KAH BN_BIH BN_ISE BN_SEB BN_DPL BN_ISW BN_RSW \
        BA_MIK_LOW BN_MIK_LOW | paste -sd '\t')

      awk -v newh="$new_headers" '
        BEGIN {{
          OFS="\t"
          print "#chr","pos","ref",newh
        }}
        {{ print }}
      ' {input.sync} > {output.sync} 2>> {log}
      """



###########################
##       READ DEPTH      ##
###########################

# rule read_depth:
#   input:
#     sync = "results/rename_sync/sync_final.sync"
#   output:
#     read_depth = "results/read_depth/read_depth.csv"
#   log:
#     "results/logs/read_depth/read_depth.log"
#   envmodules:
#     "r/4.5.0"
#   params:
#     script = "/home/eckertl/links/scratch/akpools/scripts/read_depth.R",
#     n_pops = 8,
#     sample_rate = 0.01,
#     chunk_size = 100000
#   threads: config["set-resources"]["read_depth"]["threads"]
#   shell:
#     r"""
#     mkdir -p results/read_depth results/logs/read_depth
#     Rscript {params.script} {input.sync} {params.n_pops} {params.sample_rate} {params.chunk_size} {output.read_depth} 2> {log}
#     """


###########################
##      GENOME-WIDE FST   ##
###########################

rule fst_genome:
  resources:
    mem_mb = RULE_MEM_MB["fst_genome"],
    runtime = hhmmss_to_min(RULE_TIME["fst_genome"])
  input:
    sync = "results/rename_sync/sync_final.sync",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai"
  output:
    fst_file = "results/fst_genome/fst.csv"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  log:
    "results/logs/fst_genome/fst_genome.log"
  params:
    #grenedalf_path = "/home/eckertl/software/grenedalf/bin/grenedalf",
    output_dir = "results/fst_genome",
    min_count = 2,
    min_depth = 4,
    max_depth = 200,
    pool_sizes = POOL_SIZES
  threads: config["set-resources"]["fst_genome"]["threads"]
  shell:
    r"""
    mkdir -p {params.output_dir} results/logs/fst_genome

    apptainer exec ../sif/grenedalf.sif grenedalf fst \
      --sync-path {input.sync} \
      --window-type genome \
      --window-average-policy valid-loci \
      --filter-sample-min-count {params.min_count} \
      --filter-sample-min-read-depth {params.min_depth} \
      --filter-sample-max-read-depth {params.max_depth} \
      --method unbiased-hudson \
      --pool-sizes {params.pool_sizes} \
      --threads {threads} \
      --reference-genome-fai {input.genome_fai} \
      --out-dir {params.output_dir} \
      2>> {log}
    """



###########################
##       SLIDING FST      ##
###########################

rule fst_sliding:
  resources:
    mem_mb = RULE_MEM_MB["fst_sliding"],
    runtime = hhmmss_to_min(RULE_TIME["fst_sliding"])
  input:
    sync = "results/rename_sync/sync_final.sync",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai"
  output:
    fst_file = "results/fst_sliding/fst.csv"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  log:
    "results/logs/fst_sliding/fst_sliding.log"
  params:
    #grenedalf_path = "/home/eckertl/software/grenedalf/bin/grenedalf",
    output_dir = "results/fst_sliding",
    window_width = 1000,
    min_count = 2,
    min_depth = 4,
    max_depth = 200,
    pool_sizes = POOL_SIZES
  threads: config["set-resources"]["fst_sliding"]["threads"]
  shell:
    r"""
    mkdir -p {params.output_dir} results/logs/fst_sliding

    apptainer exec ../sif/grenedalf.sif grenedalf fst \
      --sync-path {input.sync} \
      --window-type interval \
      --window-interval-width {params.window_width} \
      --window-average-policy valid-loci \
      --filter-sample-min-count {params.min_count} \
      --filter-sample-min-read-depth {params.min_depth} \
      --filter-sample-max-read-depth {params.max_depth} \
      --method unbiased-hudson \
      --pool-sizes {params.pool_sizes} \
      --threads {threads} \
      --reference-genome-fai {input.genome_fai} \
      --out-dir {params.output_dir} \
      2>> {log}
    """



###########################
##   GENOME-WIDE DIVERSITY##
###########################

rule diversity_genome:
  resources:
    mem_mb = RULE_MEM_MB["diversity_genome"],
    runtime = hhmmss_to_min(RULE_TIME["diversity_genome"])
  input:
    sync = "results/rename_sync/sync_final.sync",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai"
  output:
    diversity_file = "results/diversity_genome/diversity.csv"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  log:
    "results/logs/diversity_genome/diversity_genome.log"
  params:
    #grenedalf_path = "/home/eckertl/software/grenedalf/bin/grenedalf",
    output_dir = "results/diversity_genome",
    min_count = 2,
    min_depth = 4,
    max_depth = 200,
    pool_sizes = POOL_SIZES
  threads: config["set-resources"]["diversity_genome"]["threads"]
  shell:
    r"""
    mkdir -p {params.output_dir} results/logs/diversity_genome

    apptainer exec ../sif/grenedalf.sif grenedalf diversity \
      --sync-path {input.sync} \
      --window-type genome \
      --window-average-policy valid-loci \
      --filter-sample-min-count {params.min_count} \
      --filter-sample-min-read-depth {params.min_depth} \
      --filter-sample-max-read-depth {params.max_depth} \
      --pool-sizes {params.pool_sizes} \
      --threads {threads} \
      --reference-genome-fai {input.genome_fai} \
      --out-dir {params.output_dir} \
      2>> {log}
    """



###########################
##     SLIDING DIVERSITY  ##
###########################

rule diversity_sliding:
  resources:
    mem_mb = RULE_MEM_MB["diversity_sliding"],
    runtime = hhmmss_to_min(RULE_TIME["diversity_sliding"])
  input:
    sync = "results/rename_sync/sync_final.sync",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai"
  output:
    diversity_file = "results/diversity_sliding/diversity.csv"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  log:
    "results/logs/diversity_sliding/diversity_sliding.log"
  params:
    #grenedalf_path = "/home/eckertl/software/grenedalf/bin/grenedalf",
    output_dir = "results/diversity_sliding",
    window_width = 1000,
    min_count = 2,
    min_depth = 4,
    max_depth = 200,
    pool_sizes = POOL_SIZES
  threads: config["set-resources"]["diversity_sliding"]["threads"]
  shell:
    r"""
    mkdir -p {params.output_dir} results/logs/diversity_sliding

    apptainer exec ../sif/grenedalf.sif grenedalf diversity \
      --sync-path {input.sync} \
      --window-type interval \
      --window-interval-width {params.window_width} \
      --window-average-policy valid-loci \
      --filter-sample-min-count {params.min_count} \
      --filter-sample-min-read-depth {params.min_depth} \
      --filter-sample-max-read-depth {params.max_depth} \
      --pool-sizes {params.pool_sizes} \
      --threads {threads} \
      --reference-genome-fai {input.genome_fai} \
      --out-dir {params.output_dir} \
      2>> {log}
    """


###########################
##       RENAME SYNC     ##
###########################

rule rename_sync_nf:
  resources:
    mem_mb = RULE_MEM_MB["rename_sync_nf"],
    runtime = hhmmss_to_min(RULE_TIME["rename_sync_nf"])
  input:
    sync = "results/rm_indels/rm_indels.sync",
    bams = expand("results/rm_dupes/{pool}.bam", pool=short_names)
  output:
    sync = "results/rename_sync_nf/sync_final.sync"
  log:
    "results/logs/rename_sync_nf/rename_sync_nf.log"
  threads: config["set-resources"]["rename_sync_nf"]["threads"]
  shell:
      r"""
      mkdir -p results/rename_sync_nf results/logs/rename_sync_nf

      new_headers=$(printf "%s\n" \
        BA_PTL BA_KAH BA_BWE BA_BUN BA_RSW BA_MIK-UP BA_SEB BA_NJR BA_ISW BA_BIH \
        BA_ISE BN_NJR BN_MIK-UP BN_KAH BN_BIH BN_ISE BN_SEB BN_DPL BN_ISW BN_RSW \
        BA_MIK_LOW BN_MIK_LOW | paste -sd '\t')

      awk -v newh="$new_headers" '
        BEGIN {{
          OFS="\t"
          print "#chr","pos","ref",newh
        }}
        NR>1 {{ print }}
      ' {input.sync} > {output.sync} 2>> {log}
      """


###########################
##      GENOME-WIDE FST   ##
###########################

rule fst_genome_nf:
  resources:
    mem_mb = RULE_MEM_MB["fst_genome_nf"],
    runtime = hhmmss_to_min(RULE_TIME["fst_genome_nf"])
  input:
    sync = "results/rename_sync_nf/sync_final.sync",
    genome_fai = f"{GENOME_DIR}/{GENOME_PREFIX}.fna.fai"
  output:
    fst_file = "results/fst_genome_nf/fst.csv"
  envmodules:
    "StdEnv/2023",
    "apptainer/1.4.5"
  log:
    "results/logs/fst_genome_nf/fst_genome_nf.log"
  params:
    #grenedalf_path = "/home/eckertl/software/grenedalf/bin/grenedalf",
    output_dir = "results/fst_genome_nf",
    min_count = 2,
    min_depth = 4,
    max_depth = 200,
    pool_sizes = POOL_SIZES
  threads: config["set-resources"]["fst_genome_nf"]["threads"]
  shell:
    r"""
    mkdir -p {params.output_dir} results/logs/fst_genome_nf

    apptainer exec ../sif/grenedalf.sif grenedalf fst \
      --sync-path {input.sync} \
      --window-type genome \
      --window-average-policy valid-loci \
      --filter-sample-min-count {params.min_count} \
      --filter-sample-min-read-depth {params.min_depth} \
      --filter-sample-max-read-depth {params.max_depth} \
      --method unbiased-hudson \
      --pool-sizes {params.pool_sizes} \
      --threads {threads} \
      --reference-genome-fai {input.genome_fai} \
      --out-dir {params.output_dir} \
      2>> {log}
    """