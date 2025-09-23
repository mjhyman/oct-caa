import os
import numpy as np
import nibabel as nib
from scipy.ndimage import label
from collections import deque


def adaptive_region_growing(volume, mask, seed_thresh, intensity_tol):
    """
    Returns: Binary mask where voxels are within intensity_tol of seed points
    """
    # Only consider within mask
    vol = np.where(mask, volume, np.inf)
    seeds = np.where(vol <= seed_thresh)
    result = np.zeros(volume.shape, dtype=np.uint8)
    visited = np.zeros(volume.shape, dtype=bool)

    # 6-connectivity (for 3D)
    offsets = np.array([
        [-1, 0, 0], [1, 0, 0],
        [0, -1, 0], [0, 1, 0],
        [0, 0, -1], [0, 0, 1]
    ])

    # For every seed, grow region
    for seed_idx in zip(*seeds):
        if visited[seed_idx]:
            continue
        queue = deque()
        queue.append(seed_idx)
        seed_intensity = volume[seed_idx]
        this_region = []
        while queue:
            idx = queue.popleft()
            if visited[idx]:
                continue
            visited[idx] = True
            if (
                    abs(volume[idx] - seed_intensity) <= intensity_tol
                    and mask[idx] > 0
            ):
                result[idx] = 1
                this_region.append(idx)
                for o in offsets:
                    neighbor = tuple(np.array(idx) + o)
                    if (
                            all(0 <= n < s for n, s in zip(neighbor, volume.shape))
                            and not visited[neighbor]
                    ):
                        queue.append(neighbor)
    return result


def segment_adaptive_psoct(topdir, tissue_path, mask_path, seed_thresh, intensity_tol, min_size,
                           out_name='output_segmented.nii'):
    tissue_img = nib.load(tissue_path)
    tissue_data = tissue_img.get_fdata()
    mask_img = nib.load(mask_path)
    mask_data = mask_img.get_fdata()

    # Adaptive region growing
    print("Starting adaptive region growing...")
    rg_mask = adaptive_region_growing(tissue_data, mask_data==1, seed_thresh, intensity_tol)
    print(f"  Region growing complete. Voxels segmented: {np.sum(rg_mask)}")

    # Connected component labeling
    labeled_array, num_features = label(rg_mask)
    output = np.zeros_like(labeled_array, dtype=np.uint8)
    for n in range(1, num_features + 1):
        component = (labeled_array == n)
        if np.sum(component) >= min_size:
            output[component] = 1

    out_path = os.path.join(os.path.dirname(topdir), out_name)
    out_img = nib.Nifti1Image(output, tissue_img.affine, tissue_img.header)
    nib.save(out_img, out_path)
    print(f"Segmented mask saved to: {out_path}")


if __name__ == "__main__":
    topdir = '/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/caa17/occip'
    tissue_nifti_path = os.path.join(topdir,'caa17_occip_mus.nii')
    mask_nifti_path = os.path.join(topdir,'wm_mask_revised.nii')
    seed_threshold = 6  # intensity threshold for seeds
    intensity_tolerance = 5  # tolerance for region growing
    min_region_size = 500  # minimal size (voxels)
    segment_adaptive_psoct(topdir, tissue_nifti_path, mask_nifti_path, seed_threshold, intensity_tolerance, min_region_size)