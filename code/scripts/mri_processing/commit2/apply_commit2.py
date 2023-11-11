#!/usr/bin/python
import numpy as np
import os
import commit
from commit import trk2dictionary
import sys

#arg 1 is subj id
#arg 2 is lambda

x=sys.argv[1] #subj id

os.system( 'tck2connectome -force -nthreads 0 -assignment_radial_search 2 -out_assignments fibers_assignment.txt ' + x + '__tracking.tck ' + x + '_aparc_aseg_stages_warped.nii.gz connectome.csv' )


#Group the streamlines in bundles
if not os.path.isdir( 'bundles' ) :
    os.mkdir( 'bundles' )
os.system( 'connectome2tck -force -nthreads 0 -exclusive -files per_edge -keep_self ' + x + '__tracking.tck  fibers_assignment.txt bundles/bundle_' )

C = np.loadtxt( 'connectome.csv', delimiter=' ' ) # NB: change 'delimiter' to suits your needs
CMD = 'tckedit -force -nthreads 0'
for i in range(C.shape[0]):
    for j in range(i,C.shape[0]):
        if C[i,j] > 0 :
            CMD += ' bundles/bundle_%d-%d.tck' %(i+1,j+1)

os.system( CMD + ' '+ x + '__tracking_connecting.tck' )


#Import the tractogram

commit.core.setup()

trk2dictionary.run(
    filename_tractogram = x + '__tracking.tck',
    filename_mask       = x+ '__mask_wm.nii.gz',
    fiber_shift         = 0
)




#Initial evaluation with standard COMMIT
import amico
amico.util.fsl2scheme( 'bval', 'bvec', 'DWI.scheme' )

# load the data
mit = commit.Evaluation( '.', '.' )
mit.load_data( x + '__dwi_resampled.nii.gz', 'DWI.scheme' )


# use a forward-model with 1 Stick for the streamlines and 2 Balls for all the rest
mit.set_model( 'StickZeppelinBall' )
d_par   = 1.7E-3             # Parallel diffusivity [mm^2/s]
d_perps = [ ]                # Perpendicular diffusivity(s) [mm^2/s]
d_isos  = [ 1.7E-3, 3.0E-3 ] # Isotropic diffusivity(s) [mm^2/s]
mit.model.set( d_par, d_perps, d_isos )

mit.generate_kernels( regenerate=True )
mit.load_kernels()

# create the sparse data structures to handle the matrix A
mit.load_dictionary( 'COMMIT' )
mit.set_threads()
mit.build_operator()

# perform the fit
mit.fit( tol_fun=1e-3, max_iter=1000, verbose=True )
mit.save_results( path_suffix="_COMMIT1" )

#Retrieve the streamline contributions estimated by COMMIT for later use in COMMIT2:

x_nnls = mit.x.copy()

#Preparing the anatomical prior on bundles
C = np.loadtxt( 'connectome.csv', delimiter=' ' )
C = np.triu( C ) # be sure to get only the upper-triangular part of the matrix
group_size = C[C>0].astype(np.int32)

tmp = np.insert( np.cumsum(group_size), 0 , 0)

group_idx = np.array( [np.arange(tmp[i],tmp[i+1]) for i in range(len(tmp)-1)] )

group_w = np.empty_like( group_size, dtype=np.float64 )

for k in range(group_size.size) :
    group_w[k] = np.sqrt(group_size[k]) / ( np.linalg.norm(x_nnls[group_idx[k]]) + 1e-12 )	

#change acording to needs
reg_lambda = sys.argv[2] # change to suit your needs
#Another current limitation is the choice of the parameter ? that scales the groups penalization. From the literature we know that if the columns of the operator A are linearly independent there exists an upper (?max) an lower (?min) bounds for ? (see [24] as reference). However, this is not the case for a general tractogram, because the same (or geometrically equivalent) pathway could be shared by more than one streamline inducing redundancy inside A. In this work we set ?min = 0, which provided the results of the standard COMMIT framework, and we found the ?max empirically evaluating the loss of fibers and stopping when too many bundles were discarded resulting also in a worse fit. Future investigations to possibly set a priori a ?max, and consequently find the optimum ? are needed especially when using COMMIT on in vivo data.
prior_on_bundles = commit.solvers.init_regularisation(
    mit,
    regnorms    = [commit.solvers.group_sparsity, commit.solvers.non_negative, commit.solvers.non_negative],
    structureIC = group_idx,
    weightsIC   = group_w,
    lambdas     = [reg_lambda, 0.0, 0.0]
)
mit.fit( tol_fun=1e-3, max_iter=1000, regularisation=prior_on_bundles, verbose=True )

mit.save_results( path_suffix="_COMMIT2_" + sys.argv[2] )
