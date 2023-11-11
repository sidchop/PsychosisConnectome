#!/bin/env bash
#SBATCH --job-name=commit2
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sid.chopra@monash.edu
# SBATCH --mail-type=FAIL
# SBATCH --mail-type=BEGIN
# SBATCH --mail-type=END
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00

subj=${s}
module purge
source /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/commit2/bin/activate
module load mrtrix
unset PYTHONPATH
cd /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/commit2/${subj}
python /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/commit2/apply_commit2_lamtest.py ${subj}

