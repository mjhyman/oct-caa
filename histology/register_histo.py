import os
import glob
import SimpleITK as sitk
import numpy as np
import imageio
import matplotlib.pyplot as plt
from tkinter import filedialog, Tk

# Set this to True to save overlay images, False to only display
save_overlay = True

# Suppress root window
Tk().withdraw()

# Prompt for directory
directory = filedialog.askdirectory(title="Select Directory with OCT and Related Files")

# Find all fixed mask files
fixed_files = sorted(glob.glob(os.path.join(directory, '*_OCT_mask.tif')))

for fixed_path in fixed_files:
    prefix = os.path.basename(fixed_path).replace('_OCT_mask.tif', '')

    moving_mask_path = os.path.join(directory, f"{prefix}_5x_ds_mask.tif")
    moving_image_path = os.path.join(directory, f"{prefix}_5x_ds.tif")
    epvs_path = os.path.join(directory, f"{prefix}_5x_ds_epvs.tif")

    if not all(os.path.exists(p) for p in [moving_mask_path, moving_image_path, epvs_path]):
        print(f"Skipping {prefix}: missing one or more related files.")
        continue

    # Load images
    fixed_mask = sitk.GetImageFromArray(imageio.imread(fixed_path).astype(np.float32))
    moving_mask = sitk.GetImageFromArray(imageio.imread(moving_mask_path).astype(np.float32))
    moving_image = sitk.GetImageFromArray(imageio.imread(moving_image_path).astype(np.float32))
    epvs_image = sitk.GetImageFromArray(imageio.imread(epvs_path).astype(np.float32))

    # Registration setup
    tx = sitk.CenteredTransformInitializer(fixed_mask, moving_mask, sitk.Euler2DTransform())
    registration = sitk.ImageRegistrationMethod()
    registration.SetMetricAsMeanSquares()
    registration.SetOptimizerAsRegularStepGradientDescent(
        learningRate=2.0, minStep=1e-4, numberOfIterations=500
    )
    registration.SetInitialTransform(tx, inPlace=False)
    registration.SetInterpolator(sitk.sitkLinear)

    final_transform = registration.Execute(fixed_mask, moving_mask)

    # Resampler setup
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixed_mask)
    resampler.SetTransform(final_transform)

    # Registered moving mask
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    registered_mask = sitk.GetArrayFromImage(resampler.Execute(moving_mask)).astype(np.uint8)
    imageio.imwrite(os.path.join(directory, f"{prefix}_registered_mask.tif"), registered_mask)

    # Registered moving image
    resampler.SetInterpolator(sitk.sitkLinear)
    registered_image = sitk.GetArrayFromImage(resampler.Execute(moving_image)).astype(np.float32)
    imageio.imwrite(os.path.join(directory, f"{prefix}_registered_image.tif"), registered_image)

    # Registered EPVS image
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    registered_epvs = sitk.GetArrayFromImage(resampler.Execute(epvs_image)).astype(np.uint8)
    epvs_out_path = os.path.join(directory, f"{prefix}_registered_epvs.tif")
    imageio.imwrite(epvs_out_path, registered_epvs)

    # Overlay visualization
    fixed_np = sitk.GetArrayFromImage(fixed_mask).astype(bool)
    reg_mask_np = registered_mask.astype(bool)

    overlay = np.zeros((*fixed_np.shape, 3), dtype=np.uint8)
    overlay[..., 0] = fixed_np * 255        # Red: fixed mask
    overlay[..., 1] = reg_mask_np * 255     # Green: registered mask

    plt.figure(figsize=(6, 6))
    plt.imshow(overlay)
    plt.title(f"Overlay: {prefix}\nRed = Fixed, Green = Registered")
    plt.axis('off')
    plt.tight_layout()

    if save_overlay:
        overlay_path = os.path.join(directory, f"{prefix}_overlay.png")
        plt.savefig(overlay_path, bbox_inches='tight', dpi=150)
        print(f"Saved overlay image: {overlay_path}")

    plt.show()

    print(f"Completed registration for: {prefix}")
