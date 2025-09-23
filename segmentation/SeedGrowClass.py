import numpy as np
from collections import deque
from scipy.ndimage import gaussian_filter
from skimage.morphology import remove_small_objects


class SeedGrowingEPVS:
    def __init__(self, intensity_tolerance=10, min_size=50, smoothing_sigma=1.0):
        """
        Initialize the EPVS segmenter.

        Parameters:
        - intensity_tolerance: allowed grayscale deviation from seed intensity
        - min_size: minimum region size to keep (voxels)
        - smoothing_sigma: Gaussian blur sigma (set 0 to disable)
        """
        self.intensity_tolerance = intensity_tolerance
        self.min_size = min_size
        self.smoothing_sigma = smoothing_sigma

    def segment(self, volume, seed_mask, tissue_mask=None):
        """
        Perform 3D region-growing EPVS segmentation with an optional tissue mask.

        Parameters:
        - volume: 3D numpy array (grayscale image)
        - seed_mask: 3D binary numpy array (same shape, with manual seeds)
        - tissue_mask: 3D numpy array of same shape, with values {0,1,2}.
                       Only voxels where tissue_mask == 1 are allowed to be included.
                       If None, no tissue masking is applied.

        Returns:
        - epvs_mask: binary mask of segmented regions
        """
        if self.smoothing_sigma > 0:
            volume = gaussian_filter(volume, sigma=self.smoothing_sigma)

        X, Y, Z = volume.shape
        visited = np.zeros_like(volume, dtype=bool)
        epvs_mask = np.zeros_like(volume, dtype=bool)

        seed_coords = list(zip(*np.where(seed_mask)))

        neighbors = [(-1, 0, 0), (1, 0, 0),
                     (0, -1, 0), (0, 1, 0),
                     (0, 0, -1), (0, 0, 1)]

        for sx, sy, sz in seed_coords:
            if visited[sx, sy, sz]:
                continue

            # Check tissue_mask at seed location if provided
            if tissue_mask is not None and tissue_mask[sx, sy, sz] != 1:
                continue

            seed_intensity = volume[sx, sy, sz]
            queue = deque()
            queue.append((sx, sy, sz))

            while queue:
                x, y, z = queue.popleft()

                if not (0 <= x < X and 0 <= y < Y and 0 <= z < Z):
                    continue
                if visited[x, y, z]:
                    continue
                if tissue_mask is not None and tissue_mask[x, y, z] != 1:
                    continue

                intensity = volume[x, y, z]
                if abs(float(intensity) - float(seed_intensity)) <= self.intensity_tolerance:
                    visited[x, y, z] = True
                    epvs_mask[x, y, z] = True

                    for dx, dy, dz in neighbors:
                        xn, yn, zn = x + dx, y + dy, z + dz
                        if (0 <= xn < X) and (0 <= yn < Y) and (0 <= zn < Z):
                            if not visited[xn, yn, zn]:
                                queue.append((xn, yn, zn))

        epvs_mask = remove_small_objects(epvs_mask, min_size=self.min_size)

        return epvs_mask
