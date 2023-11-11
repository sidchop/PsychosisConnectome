# this script sets up a the need files for synb0 DISCO SDC

module load fsl

s=${subj}
dir=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi
mkdir ${dir}/syn
cp ${dir}/${s}_ses-1_dwi_denoised.nii ${dir}/syn

#extract b0 images
fslsplit ${dir}/syn/${s}_ses-1_dwi_denoised.nii ${dir}/syn/

rm ${dir}/syn/001* ${dir}/syn/002* ${dir}/syn/003* ${dir}/syn/004* ${dir}/syn/005* ${dir}/syn/006* ${dir}/syn/${s}_ses-1_dwi_denoised.nii -rf

fslmerge -t ${dir}/syn/b0 ${dir}/syn/000*
rm ${dir}/syn/000* -rf

#create  mean b0
fslmaths ${dir}/syn/b0.nii.gz -Tmean ${dir}/syn/b0.nii.gz

#copy acqprams from sub-001
cp /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/sub-001/ses-1/dwi/syn/acqparams.txt ${dir}/syn

#copy and rename T1 
cp /home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/anat/${s}_ses-1_T1w.nii ${dir}/syn 
mv ${dir}/syn/${s}_ses-1_T1w.nii ${dir}/syn/T1.nii 
gzip ${dir}/syn/T1.nii


#done
echo ${s}
