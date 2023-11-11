#Denoise - Step 1 of pipeline

module load mrtrix
dwidenoise /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi.nii /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi_denoised.nii
