# -*- coding: utf-8 -*-
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap

# File paths
cwrf_file = '/share/home/dq135/CRESM/CASES/ICBC_CN01_2000/wrfinput_d01'
colm_file = '/share/home/dq135/draw/scripts/land_sea_distribution/unstructured_cwrf_CN01_hist_2000-05_remap.nc'

# Custom colormap for only 2 and 3
colors = ['green', 'red']  # 2: green, 3: red
cmap = ListedColormap(colors)
bounds = [2, 2.5, 3]  # Define boundaries to map 2 and 3 to colors
norm = plt.Normalize(2, 3)

# Read CWRF data
with nc.Dataset(cwrf_file, 'r') as ds_cwrf:
    landmask_cwrf = ds_cwrf.variables['LANDMASK'][0, :, :]  # 1: land, 0: water
    lat_cwrf = ds_cwrf.variables['XLAT'][0, :, :]  # 2D latitude
    lon_cwrf = ds_cwrf.variables['XLONG'][0, :, :]  # 2D longitude

# Read CoLM data
with nc.Dataset(colm_file, 'r') as ds_colm:
    f_t_soisno = ds_colm.variables['f_t_soisno'][:]  # Read all data
    print("f_t_soisno shape:", f_t_soisno.shape)  # Debug: print initial shape
    f_t_soisno_slice = f_t_soisno[0, :, :, :]  # Take first time step
    print("f_t_soisno[0, :, :, :] shape:", f_t_soisno_slice.shape)  # Debug: print sliced shape
    # Average over the last axis dynamically
    try:
        f_t_soisno_2d = np.nanmean(f_t_soisno_slice, axis=-1)  # Average over the last dimension
    except ValueError as e:
        print(f"Error in nanmean: {e}. Using raw slice instead.")
        f_t_soisno_2d = f_t_soisno_slice
    # Define landmask: 1 where f_t_soisno has valid values, 0 where NaN or -9999
    landmask_colm = np.where(np.isnan(f_t_soisno_2d) | (f_t_soisno_2d == -9999), 0, 1)

# Verify masks
print("CWRF LANDMASK unique:", np.unique(landmask_cwrf))
print("CoLM LANDMASK unique:", np.unique(landmask_colm))

# Ensure grid sizes match
assert landmask_cwrf.shape == landmask_colm.shape, "The grid sizes of the two files do not match!"

# Create difference mask (only 2 and 3)
diff_mask = np.full(landmask_cwrf.shape, np.nan, dtype=float)
diff_mask[(landmask_cwrf == 1) & (landmask_colm == 0)] = 2  # CWRF land, CoLM ocean
diff_mask[(landmask_cwrf == 0) & (landmask_colm == 1)] = 3  # CWRF ocean, CoLM land

# Print mask distribution for verification
print("DIFF_MASK unique values:", np.unique(diff_mask[~np.isnan(diff_mask)], return_counts=True))

# Visualization
plt.figure(figsize=(12, 8))
proj = ccrs.LambertConformal(
    central_longitude=110.0,  # CEN_LON
    central_latitude=35.1778,  # CEN_LAT
    standard_parallels=(30.0, 60.0)  # TRUELAT1, TRUELAT2
)
ax = plt.axes(projection=proj)

# Add China boundary (offline mode, replace with your local Shapefile path)
try:
    import cartopy.io.shapereader as shpreader
    china_shp = shpreader.Reader('/share/home/dq135/draw/shapefile_China/ne_110m_coastline.shp')  # Replace with actual path
    ax.add_geometries(china_shp.geometries(), ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.5)
except (ImportError, ValueError):
    print("Warning: Could not load China Shapefile. Using empty background.")

img = ax.pcolormesh(lon_cwrf, lat_cwrf, diff_mask, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())

# Add colorbar with only 2 and 3 labels
cbar = plt.colorbar(img, ticks=[2, 3])
cbar.ax.set_yticklabels(['CWRF Land, CoLM Ocean', 'CWRF Ocean, CoLM Land'])

# Set title and labels
plt.title('CWRF and Hawkins CoLM Land-Sea Distribution Difference')
plt.savefig('sea_land_diff.pdf', bbox_inches='tight', format='pdf')
plt.close()

print("Difference mask created and visualized. Image saved as 'sea_land_diff.pdf'.")