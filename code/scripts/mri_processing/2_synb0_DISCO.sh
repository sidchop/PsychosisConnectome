#!/bin/sh
#SBATCH --job-name=STAGES_fmriprep
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
# SBATCH --mail-user=sid.chopra@monash.edu
# SBATCH --mail-type=BEGIN
# SBATCH --mail-type=FAIL
# SBATCH --mail-type=END
#SBATCH --mem-per-cpu=8000
#SBATCH --qos=normal


module load singularity

s=${subj}

singularity run -e \
-B /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/syn/:/INPUTS \
-B /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/:/OUTPUTS \
-B /home/scho0011/kg98/Sid/license.txt:/extra/freesurfer/license.txt \
/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/synb0_latest.sif \

