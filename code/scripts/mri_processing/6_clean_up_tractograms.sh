#remove random streamlines outside brain
module purge
module load mrtrix
for s in `cat temp_list.txt` ; do tckedit /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${s}/all.tck /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/commit2/${s}/all_cropped.tck -mask /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow_pve/results/${s}/Register_T1/${s}__t1_mask_warped.nii.gz -force ; echo ${s} ; done


