mkdir /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${subj} ; cp /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/tractoflow_root/${subj}/bvec /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/tractoflow_root/${subj}/bval /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${subj} ; cp /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_pve/results/${subj}/Tracking/all.trk /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_pve/results/${subj}/Resample_DWI/${subj}__dwi_resampled.nii.gz /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_pve/results/${subj}/FODF_Metrics/${subj}__peaks.nii.gz /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_pve/results/${subj}/Register_T1_Maps/${subj}__mask_wm_warped.nii.gz /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${subj} ; cp /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/atlas/included_rois_only/${subj}_schaefer300n7_aseg_to_dwispace_gm_rois.nii.gz /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${subj}

module purge
module load nibabel/2.3.3

python trk2tck.py /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${subj}/all.trk -f
if [ -f /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${subj}/all.tck ] ; then rm /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${subj}/all.trk -f ; fi

echo ${subj}
