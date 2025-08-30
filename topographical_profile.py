# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 11:01:10 2025

@author: vigno
"""

# ===============================
# Import needed libraries
# ===============================

import rasterio
import numpy as np
import matplotlib.pyplot as plt
from rasterio.transform import rowcol
from math import cos, radians
from scipy.spatial.distance import euclidean

# ===============================
# Definition of functions
# ===============================
def bresenham_line(x0, y0, x1, y1):
    """
    Draw the line between two points.
    """
    points = []
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    x, y = x0, y0
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx - dy
    while True:
        points.append((x, y))
        if x == x1 and y == y1:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x += sx
        if e2 < dx:
            err += dx
            y += sy
    return points


def get_profile(raster_path, lon1, lat1, lon2, lat2, plot_map=True):
    """
    Compute a topographical profile between two points and shows (1) the map with the profile and (2) the topographical profile.
    
    Parameters
    ----------
    raster_path : str
        Path
    lon1, lat1, lon2, lat2 : float
        Geographical or projected coordinates
    plot_map : bool
        If True, plot the map with the profile
    
    Output
    --------
    distances, elevations
    """
    with rasterio.open(raster_path) as src:
        raster = src.read(1)
        transform = src.transform
        crs = src.crs
        nrows, ncols = raster.shape

        row1, col1 = rowcol(transform, lon1, lat1)
        row2, col2 = rowcol(transform, lon2, lat2)

        # Récupération des bornes du raster
        bounds = src.bounds  # left, bottom, right, top

        if not (0 <= row1 < nrows and 0 <= col1 < ncols):
            raise ValueError(
                f"Point A ({lon1}, {lat1}) is not in the raster extent!\n"
                f"Raster bounds:\n"
                f"  Longitude/X: {bounds.left:.6f} to {bounds.right:.6f}\n"
                f"  Latitude/Y : {bounds.bottom:.6f} to {bounds.top:.6f}"
                )

        if not (0 <= row2 < nrows and 0 <= col2 < ncols):
            raise ValueError(
                f"Point B ({lon2}, {lat2}) is not in the raster extent!\n"
                f"Raster bounds:\n"
                f"  Longitude/X: {bounds.left:.6f} to {bounds.right:.6f}\n"
                f"  Latitude/Y : {bounds.bottom:.6f} to {bounds.top:.6f}"
                )


        # Pixels of the profile, defined in the bresenham_line function
        pixels = bresenham_line(col1, row1, col2, row2)
        valid_pixels = [(c, r) for c, r in pixels if 0 <= r < nrows and 0 <= c < ncols]

        elevations = [raster[r, c] for c, r in valid_pixels]

        # Distances computations
        distances = [0.0]
        if crs and crs.is_geographic:
            lat_mean = (lat1 + lat2) / 2.0
            deg2km_lat = 111.32                             # for Earth
            deg2km_lon = 111.32 * cos(radians(lat_mean))    # for Earth
            for i in range(1, len(valid_pixels)):
                lon_px, lat_px = src.xy(valid_pixels[i][1], valid_pixels[i][0])
                lon_prev, lat_prev = src.xy(valid_pixels[i-1][1], valid_pixels[i-1][0])
                dx = (lon_px - lon_prev) * deg2km_lon
                dy = (lat_px - lat_prev) * deg2km_lat
                distances.append(distances[-1] + np.sqrt(dx**2 + dy**2))
            xunit = "Distance (km)"
        else:
            for i in range(1, len(valid_pixels)):
                d_pix = euclidean(valid_pixels[i], valid_pixels[i - 1])
                d_m = d_pix * abs(transform.a)
                distances.append(distances[-1] + d_m)
            xunit = "Distance (m)"

        # --- map with profile ---
        fig, ax = plt.subplots(1, 2 if plot_map else 1, figsize=(12 if plot_map else 8, 4))

        if plot_map:
            ax_map, ax_prof = ax
            extent = rasterio.plot.plotting_extent(src)
            im = ax_map.imshow(raster, extent=extent, cmap="terrain", origin="upper")
            cbar = plt.colorbar(im, ax=ax_map, orientation="horizontal", fraction=0.05, pad=0.15)
            cbar.set_label("Elevation") 

            xs, ys = zip(*[src.xy(r, c) for c, r in valid_pixels])
            ax_map.plot(xs, ys, "r--", linewidth=2, label="Profile")
            ax_map.scatter([lon1, lon2], [lat1, lat2], c="k", marker=".", zorder=5)
            ax_map.set_title("Map")
            ax_map.legend()
        else:
            ax_prof = ax

        ax_prof.plot(distances, elevations, "-k")
        ax_prof.set_xlabel(xunit)
        ax_prof.set_ylabel("Elevation")
        ax_prof.set_title("Topographic profile")
        ax_prof.grid(True)

        plt.tight_layout()
        plt.show()

        return distances, elevations

# ===============================
if __name__ == "__main__":
    raster_file = "your_raster.tif"  # name/directory of your raster (.tif, .tiff, .png. jpg, .nc)

    # Two points that are inside the raster extent
    lon1, lat1 = 0, -60
    lon2, lat2 = 170, -60

    try:
        dist, elev = get_profile(raster_file, lon1, lat1, lon2, lat2)

    except ValueError as e:
        print("Error:", e)
