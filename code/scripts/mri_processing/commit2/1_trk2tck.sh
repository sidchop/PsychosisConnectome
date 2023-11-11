module load purge
module load nibable

python trk2tck.py TractoFlow/results/${subj}/Tracking/${subj}__tracking.trk -f

module load mrtrix 
tck2connectome -force -nthreads 0 -assignment_radial_search 2 -out_assignments /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/commit2/${subj}/fibers_assignment.txt TractoFlow/results/${subj}/Tracking/${subj}__tracking.tck /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/aseg_aparc_atlas/${subj}_aparc_aseg_stages_warped.nii.gz /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/commit2/${subj}/connectome.csv

