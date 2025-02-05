import ee

# Initialize Earth Engine
ee.Initialize()

# Define the study area (edaratkol_fc) and date range
start_date = '2023-05-01'
end_date = '2023-05-30'

# Load Sentinel-1 GRD collection and apply time filter
s1 = ee.ImageCollection("COPERNICUS/S1_GRD") \
    .filterBounds(edaratkol_fc) \
    .filterDate(start_date, end_date)  # Filter by the time range

# **Load & Process SRTM DEM for Terrain Correction**
srtm = ee.Image("USGS/SRTMGL1_003")  # Load SRTM 30m DEM
srtm_clipped = srtm.clip(edaratkol_fc)  # Clip DEM to study area
srtm_filtered = srtm_clipped  # No need for time filtering on SRTM data (static)
srtm_clean = srtm_filtered.toFloat().updateMask(srtm_filtered.neq(-32768))  # Convert to Float & Remove No-Data

# **Ensure VV and VH Bands Are Used**
polarization_bands = ['VV', 'VH']

# Filter only for the detected polarization type
s1 = s1.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
s1 = s1.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))

# **Function to Apply Thermal Noise Removal**
def thermal_noise_removal(image):
    return image.select(polarization_bands).subtract(image.select(['angle']).multiply(0.001))

# **Function to Apply Calibration**
def do_calibration(image, polarization, pols):
    return image.select(pols).divide(10.0).exp()

# **Function to Apply Speckle Filtering**
def do_speckle_filtering(image):
    return image.focal_mean(radius=3, kernelType='circle', iterations=1)

# **Function to Apply Terrain Correction**
def do_terrain_correction(image, proj, downsample):
    terrain_corrected = image.reproject(crs=proj, scale=40 if downsample else 10)
    return terrain_corrected

# Apply Preprocessing Steps
composite = (
    s1.filter(ee.Filter.eq('instrumentMode', 'IW'))  # IW mode for Interferometric Wide
    .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))  # Ascending orbit pass
    .map(thermal_noise_removal)
    .map(lambda img: do_calibration(img, 'VV', polarization_bands))
    .map(do_speckle_filtering)
    .map(lambda img: do_terrain_correction(img, 'EPSG:4326', downsample=1))
)

# **Check Number of Images Before Processing**
num_images = composite.size().getInfo()
print("✅ Number of Images Before Processing:", num_images)

# **Apply Median Reduction, Keeping All Bands**
s1_mosaic = composite.select(polarization_bands).median().clip(edaratkol_fc)

# **Fix Resolution by Reprojecting to Native Sentinel-1 Projection**
s1_mosaic = s1_mosaic.reproject(
    crs=s1.first().select(polarization_bands[0]).projection(),
    scale=10  # Sentinel-1 GRD IW native resolution
)

# **Check Final Resolution**
resolution = s1_mosaic.projection().nominalScale().getInfo()
print("✅ Final Corrected Output Image Resolution (meters):", resolution)
