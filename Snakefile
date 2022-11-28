#! /usr/bin/py
import pandas as pd 
import os

#trim -> fastqc -> mapping -> sort bam -> index -> bigWig

configfile: "config.json"
localrules: all, mkdir 

df = pd.read_csv(config["meta_file"], sep='\t', header=0, index_col=0)
sample_ids = list(df.index)
print(df.index)

def get_gz(sample_id):
    dir=config["raw_fastq_gz_dir"]
    print(sample_id) 
    print(dir)
    print (os.path.join(dir,df.loc[sample_id]["ForwardFastqGZ"]))
    return os.path.join(dir,df.loc[sample_id]["ForwardFastqGZ"])
    
def get_forward_primer(sample_id):
    return df.loc[sample_id]["Adapter_1"]

rule all:
    input:expand("{dir}/{sample_id}.cpm.norm.bw",dir=config["dir_names"]["bigwigs_dir"], sample_id=sample_ids)
    run:
        for sample in sample_ids:
            print("Wrapping up pipeline")

rule mkdir:
    output: touch(config["file_names"]["mkdir_done"])
    params: dirs = list(config["dir_names"].values())
    resources: time_min=10, mem_mb=2000, cpus=1
    shell: "mkdir -p {params.dirs}"

rule trim: 
    input:
        rules.mkdir.output,
        all_read1 = lambda wildcards: get_gz(wildcards.sample_id)
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        trimmed_read1 = config["dir_names"]["trimmed_dir"]+ "/{sample_id}.trimmed.R1.fastq.gz",
        trimmed_stats = config["dir_names"]["trimmed_dir"]+ "/{sample_id}.trimmed.stats"
    version: config["tool_version"]["cutadapt"]
    params:
        adapter1=lambda wildcards: get_forward_primer(wildcards.sample_id)
    shell: "cutadapt -m 15 -a {params.adapter1} -n 1 -o {output.trimmed_read1} {input.all_read1} >{output.trimmed_stats}"

rule fastqc:
    input:
        read1= rules.trim.output.trimmed_read1
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        fastqc_file1 = config["dir_names"]["fastqc_dir"]+ "/{sample_id}.trimmed.R1_fastqc.html",
    shell: "fastqc {input.read1} -q -o outputs/fastqc/" 

rule map:
    input:
        p1 = rules.trim.output.trimmed_read1,
        p2 = rules.fastqc.output.fastqc_file1
    resources: time_min=360, mem_mb=20000, cpus=6
    output: 
        mapped_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam",
        stats = config["dir_names"]["mapped_dir"]+"/{sample_id}.stats",
    params:
        threads = config["params"]["bowtie2"]["threads"],
        map_all = config["params"]["bowtie2"]["all"],
        reference = config["params"]["bowtie2"]["bowtie2_reference"]
   # shell:"bowtie2 --threads {params.threads} -1 {input.p1} -x {params.reference} 2> {output.stats} | samtools sort -T {output.mapped_bam_file}.tmp -O bam -o {output.mapped_bam_file}"
    shell:
        """
        bowtie2 --threads {params.threads} {input.p1} -x {params.reference} 2> {output.stats} | samtools sort -T {output.mapped_bam_file}.tmp -O bam -o {output.mapped_bam_file}
        """

rule sort_bam:
    input:
        sorted_bam = rules.map.output.mapped_bam_file
    output:
        sorted_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.sorted.bam"
    resources: time_min=360, mem_mb=20000, cpus=6
    shell: "samtools view -hb {input.sorted_bam} | samtools sort -T {input.sorted_bam}.tmp -o {output.sorted_bam_file}" 

rule index_bam:
    input:
        index_bam = rules.sort_bam.output.sorted_bam_file
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        index_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.sorted.bam.bai"
    shell: "samtools index {input.index_bam} {output.index_bam_file}"

rule normalized_bigwig:
    input:
        sorted_bams = rules.sort_bam.output.sorted_bam_file,
        index_bams = rules.index_bam.output.index_bam_file
    resources: time_min=360, mem_mb=20000, cpus=6
    output:
        normalized_bigwig_file = config["dir_names"]["bigwigs_dir"] + "/{sample_id}.cpm.norm.bw"
    params:
        cores = config["params"]["deeptools"]["cores"],
        size = config["params"]["deeptools"]["size"]
    shell:
        """
        bamCoverage --bam {input.sorted_bams} -o {output.normalized_bigwig_file} \
            --binSize 10 \
            --normalizeUsing CPM \
            --effectiveGenomeSize {params.size} \
            --numberOfProcessors {params.cores}
        """
