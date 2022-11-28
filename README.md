# ChIP-Seq
ChIP-Seq snakemake for single-end sequences

Currently Loaded Modules:

`module load gbc-samtools/1.12 gbc-deeptools/3.4.3 gbc-bowtie2/2.4.2 gbc-cutadapt/1.16 gbc-bedtools/2.29.1 python/py37-anaconda-2019.10 snakemake/5.7.1-py37`

Step-by-step of install and analysis

1. git clone this repository

`git clone https://github.com/achisha21/ChIP-Seq.git`

2. Activate the python anaconda environment

`conda activate snakemake`

3. Edit the config.json file and cluster.json files

4. Ensure meta-data table contains all of the necessairy fields

** NOTE EXACT HEADERS HAVE TO BE ENFORCED or key errors will be thrown during processing**

5. Launch jobs
The use of --latency-wait allows for SLURM to catch up writing the files and posting the file handles so Snakemake can see them.

`snakemake --latency-wait 120 -p -j 100 --profile slurm`

6. Pipeline should result in bigwig files

#Using gnu parallel for running MACS2 for peak calling 

Load the following modules

`module load gnu-parallel MACS/2.2.7.1`

`for bam in $(ls *.sorted.bam); do echo $bam; done`

`for bam in $(ls *.sorted.bam); do echo "macs2 callpeak -t $bam -n $bam.tmp -g hs --outdir macs2 "; done >> cmds.txt`

`parallel -j 8 < example.cmds.txt`

Run this from an interactive job (Not on vortex)

`parallel -j 8 < cmds.txt`


