#!/bin/bash
#SBATCH --time=14-0:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --job-name bt2
#SBATCH --mail-type=ARRAY_TASKS,END
#SBATCH --mail-user=
#SBATCH --output=%x-%j.SLURMout

module load Bowtie2
module load samtools/1.12-rhel8

# SLURM Array Config
config=/work/[netID]/knockout/array_ids.txt
i=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
echo -e "${i}_${SLURM_ARRAY_TASK_ID} ${SLURM_JOB_NAME}\n"

# aedes aegypti Bowtie2-indexed reference genome
REF_GENOME="/hpc/group/[lab]/ref/genomes/a_aegypti/GCF_002204515.2_AaegL5.0-Bowtie2-build"

cd -- /work/[netID]/knockout/$i

bowtie2 -x $REF_GENOME -q -1 $i\_R1.fastq.gz -2 $i\_R2.fastq.gz | samtools view -bS - | samtools sort - -o $i\_sorted.bam

# index sorted bam file and generate a bai file
samtools index -b $i\_sorted.bam
echo "Done!"
