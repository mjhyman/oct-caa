import nibabel as nib
import numpy as np
from SeedGrowClass import SeedGrowingEPVS

# Load volume and seeds
volume_nifti = nib.load('/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/'
                        'caa17/occip/caa17_occip_mus.nii')
seeds_nifti = nib.load('/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/'
                        'caa17/occip/epvs.nii')
tissue_mask = nib.load('/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/'
                        'caa17/occip/wm_mask_revised.nii').get_fdata().astype(np.int8)  # expected values: 0,1,2

# Convert to data types
volume = volume_nifti.get_fdata().astype(np.float32)
seed_mask = seeds_nifti.get_fdata().astype(bool)

# Create segmenter
segmenter = SeedGrowingEPVS(intensity_tolerance=10,min_size=50,smoothing_sigma=2.0)

# Run segmentation
epvs_mask = segmenter.segment(volume, seed_mask, tissue_mask=tissue_mask)

# Save result
nib.save(nib.Nifti1Image(epvs_mask.astype(np.uint8), affine=np.eye(4)), '/autofs/cluster/octdata3/users/'
                                                                        'mjhyman/oct_caa_analyses/optical_properties/'
                                                                        'caa17/occip/epvs_seed_grow.nii')