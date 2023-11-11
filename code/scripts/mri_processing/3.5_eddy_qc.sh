#run quad

module load fsl
s=${SUBJ}

eddy_quad /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi_denoised_eddy -idx /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/index.txt -par /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/syn/acqparams.txt -m /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi_bet_mask.nii.gz -b /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi/${s}_ses-1_dwi.bval -s /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/my_slspec.txt -v

