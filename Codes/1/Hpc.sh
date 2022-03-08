#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=Hpc.sh
#SBATCH --output=Hpc.sh.o%j
#SBATCH --error=Hpc.sh.e%j
#SBATCH --time=36:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c.liu@oist.jp

module load R/4.0.4
module load python/3.10.2

module load bioinfo-ugrp-modules
module load DebianMed/11.2

module load sra-toolkit/2.10.9+dfsg-2
module load fastqc/0.11.9
module load Trimmomatic/0.39
module load hisat2/2.2.1-2+b3
module load samtools/1.12
module load stringtie/2.1.4+ds-4

cd /home/c/c-liu/DEG

Rscript Run.R metadata.txt . ReferenceGenome/NemVec ReferenceGenome/NemVec.gff 32