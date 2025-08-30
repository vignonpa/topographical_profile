# Topographic Profile Extraction

## Description

This Python script allows you to extract **topographic profiles** from a georeferenced raster (e.g., DEM, GeoTIFF) between two points. It works on Earth (by default) but also on other planetary bodies when changing the deg_to_km values. It provides:  

- A plot of **elevation vs. distance** along the profile.  
- An optional **map view** of the raster with the profile line and start/end points, and an associated colorbar.  
- Automatic **distance conversion to kilometers** if the raster is in geographic coordinates (longitude/latitude).  
- **Boundary checks** to ensure points are inside the raster extent.  

This tool is useful for **geosciences, planetary science, and earth/planetary data analysis**.

---

## Requirements

- Python 3.8+  
- [rasterio](https://rasterio.readthedocs.io/en/latest/)  
- [numpy](https://numpy.org/)  
- [matplotlib](https://matplotlib.org/)  
- [scipy](https://www.scipy.org/)  

Install dependencies via pip:

```bash
pip install rasterio numpy matplotlib scipy
