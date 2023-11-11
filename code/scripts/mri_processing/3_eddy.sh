#!/bin/bash
#SBATCH --job-name=Eddy
#SBATCH --account=kg98
#SBATCH --time=56:00:00
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --partition=m3h

module load cuda/8.0
module load fsl

s=${subj}

eddy_cuda8.0 \
--imain=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi_denoised.nii \
--mask=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi_bet_mask.nii.gz \
--acqp=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/syn/acqparams.txt \
--index=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/index.txt \
--bvecs=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi.bvec \
--bvals=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi.bval \
--topup=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/topup \
--out=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi_denoised_eddy \
--cnr_maps \
--repol \
--mporder=16 \
--s2v_niter=10 \
--s2v_lambda=5 \
--s2v_interp=trilinear \
--slspec=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/my_slspec.txt \
--verbose

