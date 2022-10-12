#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name index
#SBATCH --output=%x-%j.SLURMout

#Change to current working directory
cd ${PBS_O_WORKDIR}

#Add conda environment to the path
export PATH="${HOME}/miniconda3/envs/plb812/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/miniconda3/envs/plb812/lib:${LD_LIBRARY_PATH}"

hisat2-build -p 16 --exon Sorghum_bicolor_NCBIv3_genomic.exon --ss Sorghum_bicolor_NCBIv3_genomic.ss Sorghum_bicolor_NCBIv3_genomic.fna  sorghum_index
