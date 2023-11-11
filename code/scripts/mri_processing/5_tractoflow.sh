#!/bin/sh

#SBATCH --job-name=STAGES_tractoflow
#SBATCH --account=kg98
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=48:00:00

# export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)

module load singularity
sh /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_cont/nextflow run /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_cont/tractoflow-2.1.1/main.nf --root /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/tractoflow_root/cont/ --dti_shells "0 3000" --fodf_shells "0 3000" --run_dwi_denoising "false" --run_topup "false" --algo "det" --processes_brain_extraction_t1 8 --processes_denoise_t1 8 --processes_fodf 8 --processes_registration 8 --run_eddy "false" -with-singularity /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_cont/tractoflow_2.1.1_650f776_2020-07-15.img 
