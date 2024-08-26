import os

import numpy as np
import pyelastix
import matplotlib.pyplot as plt

# os.environ["ELASTIX_PATH"] = "/home/nico/snap/elastix-5.0.1-linux/bin/elastix"
# os.environ["LD_LIBRARY_PATH"] = "/home/nico/snap/elastix-5.0.1-linux/lib"

os.environ["ELASTIX_PATH"] = "elastix/elastix-5.1.0/bin/"
os.environ["LD_LIBRARY_PATH"] = "elastix/elastix-5.1.0/lib/"

class ElastixRegistration:

    def __init__(self, params: dict) -> None:
        self.params = pyelastix.get_default_params()
        self.params.MaximumNumberOfIterations = params['max_num_it']
        self.params.NumberOfResolutions = params['nr_res']
        self.params.Metric = params['metric']
        self.params.NumberOfHistogramBins = params['nr_histogram_bins']
        self.params.NumberOfSpatialSamples = params['nr_spatial_samples']
        self.params.FinalGridSpacingInPhysicalUnits = params['bspline_grid_spacing']
        self.params.Transform = params["transform"]


    def register_images(self, fixed: np.ndarray, moving: np.ndarray, plot_vectorfield=False, return_vectorfield=False) -> np.ndarray:
        moving = np.ascontiguousarray(moving, dtype=np.float32)
        fixed = np.ascontiguousarray(fixed, dtype=np.float32)

        deformed, field = pyelastix.register(moving, fixed, self.params, verbose=0)
        deformed = np.clip(deformed, 0, 255).astype(np.uint8)

        if plot_vectorfield:
            v, u = field
            self.plot_vectorfield(u, v, fixed)

        return (deformed, field) if return_vectorfield else deformed


    def plot_vectorfield(self, u, v, fixed):
        norm = np.sqrt(u ** 2 + v ** 2)
        nvec = 40  # Number of vectors to be displayed along each image dimension
        nl, nc = fixed.shape
        step = max(nl // nvec, nc // nvec)

        y, x = np.mgrid[:nl:step, :nc:step]
        u_ = u[::step, ::step]
        v_ = v[::step, ::step]

        fig, axs = plt.subplots()
        imshow = axs.imshow(norm)
        plt.quiver(x, y, u_, v_, color='r', units='dots', angles='xy', scale_units='xy', lw=3)
        plt.colorbar(imshow, ax=axs, fraction=0.046, pad=0.04, label="displacement magnitude")
        axs.set_title(f"Vector field magnitude and direction for {self.__class__.__name__}")
