cd /home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet

module purge
module load singularity
singularity shell ../set_1v0/set_1v0.img


for s in `cat subjlist.txt` ; do    \

t1=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}/final/*_t1_final.mnc 
xfm_transfo=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}/transforms/linear/*t1_tal.xfm 
pve_csf=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}/classify/*pve_exactcsf.mnc 
pve_gm=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}/classify/*pve_exactgm.mnc 
pve_wm=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}/classify/*pve_exactwm.mnc 
pve_sc=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}/classify/*pve_exactsc.mnc 
brain_mask=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}/mask/*brain_mask.mnc 

mincresample -trilinear -tfm_input_sampling -invert_transformation -transformation $xfm_transfo $t1 -clobber ${s}__t1.mnc ;
mnc2nii ${s}__t1.mnc ${s}__temp.nii ;
scil_resample_volume.py ${s}__temp.nii ${s}__t1.nii.gz --ref ${s}__temp.nii --enforce_dimensions ;

mincresample -trilinear -like ${s}__t1.mnc -tfm_input_sampling -invert_transformation -transformation $xfm_transfo $pve_csf -clobber ${s}__temp.mnc ;
mnc2nii ${s}__temp.mnc ${s}__temp.nii ;
scil_resample_volume.py ${s}__temp.nii ${s}__pve_csf.nii.gz --ref ${s}__t1.nii.gz --enforce_dimensions ;

mincresample -trilinear -like ${s}__t1.mnc -tfm_input_sampling -invert_transformation -transformation $xfm_transfo $pve_gm -clobber ${s}__temp.mnc ;
mnc2nii ${s}__temp.mnc ${s}__temp.nii ;
scil_resample_volume.py ${s}__temp.nii ${s}__pve_gm.nii.gz --ref ${s}__t1.nii.gz --enforce_dimensions ;

mincresample -trilinear -like ${s}__t1.mnc -tfm_input_sampling -invert_transformation -transformation $xfm_transfo $pve_sc -clobber ${s}__temp.mnc ;
mnc2nii ${s}__temp.mnc ${s}__temp.nii ;
scil_resample_volume.py ${s}__temp.nii ${s}__pve_sc.nii.gz --ref ${s}__t1.nii.gz --enforce_dimensions ;

mincresample -trilinear -like ${s}__t1.mnc -tfm_input_sampling -invert_transformation -transformation $xfm_transfo $pve_wm -clobber ${s}__temp.mnc ;
mnc2nii ${s}__temp.mnc ${s}__temp.nii ;
scil_resample_volume.py ${s}__temp.nii ${s}__pve_wm.nii.gz --ref ${s}__t1.nii.gz --enforce_dimensions ;

mincresample -nearest_neighbour -like ${s}__t1.mnc -tfm_input_sampling -invert_transformation -transformation $xfm_transfo $brain_mask -clobber ${s}__temp.mnc \
mnc2nii ${s}__temp.mnc ${s}__temp.nii ;
scil_resample_volume.py ${s}__temp.nii ${s}__brain_mask.nii.gz --ref ${s}__t1.nii.gz --enforce_dimensions --interp 'nn' ;

rm ${s}__temp.mnc ${s}__temp.nii ${s}__t1.mnc ; done

exit #exit singularity shell


#move to tactoflow folder 

for s in `cat subjlist.txt` ; do \
tracto_input=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/tractoflow_root/${s}
t1=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}__t1.nii.gz
pve_csf=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}__pve_csf.nii.gz
pve_gm=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}__pve_gm.nii.gz
pve_wm=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}__pve_wm.nii.gz
pve_sc=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}__pve_sc.nii.gz
brain_mask=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/civet/${s}__brain_mask.nii.gz

cp $t1 $pve_csf $pve_gm $pve_wm $pve_sc $brain_mask temp ;
mv temp/${s}__t1.nii.gz temp/t1.nii.gz ;
mv temp/${s}__pve_csf.nii.gz temp/pve_csf.nii.gz ;
mv temp/${s}__pve_gm.nii.gz temp/pve_gm.nii.gz ;
mv temp/${s}__pve_wm.nii.gz temp/pve_wm.nii.gz ;
mv temp/${s}__pve_sc.nii.gz temp/pve_sc.nii.gz ;
mv temp/${s}__brain_mask.nii.gz temp/brain_mask.nii.gz ;

mv temp/* $tracto_input ; echo ${s} ; done



