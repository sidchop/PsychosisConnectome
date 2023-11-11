#!/bin/bash

s=${SUBJ}
module load fsl
dir=/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/${s}/ses-1/dwi
bet ${dir}/${s}_ses-1_dwi.nii ${dir}/${s}_ses-1_dwi_bet.nii -m -R -f 0.16
