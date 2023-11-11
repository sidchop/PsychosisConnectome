s=${SUBJ}

module load freesurfer

mri_convert  /home/scho0011/kg98/Sid/STAGES/freesurfer/${s}.long.${s}_base/mri/aparc+aseg.mgz /home/scho0011/kg98/Sid/STAGES/freesurfer/${s}.long.${s}_base/mri/aparc+aseg.nii.gz

mri_convert  /home/scho0011/kg98/Sid/STAGES/freesurfer/${s}.long.${s}_base/mri/wmparc.mgz /home/scho0011/kg98/Sid/STAGES/freesurfer/${s}.long.${s}_base/mri/wmparc.nii.gz

mv /home/scho0011/kg98/Sid/STAGES/freesurfer/${s}.long.${s}_base/mri/aparc+aseg.nii.gz /home/scho0011/kg98_scratch/Sid/STAGES_dti/tractoflow_abs_root/sub-${s}

mv /home/scho0011/kg98/Sid/STAGES/freesurfer/${s}.long.${s}_base/mri/wmparc.nii.gz /home/scho0011/kg98_scratch/Sid/STAGES_dti/tractoflow_abs_root/sub-${s}
