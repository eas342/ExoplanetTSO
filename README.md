# ExoplanetTSO
# Exoplanet Photometric Analysis for Time Series Observations
---

Many ground and space telescopes missions contain the capacity to perform photometric time series observations.  
In order to process that data, it is important to identify the host star, accurately estimate it's PSF (point spread function) 
center -- for all images in the time series.

## JWST - NIRCam - Level 3 TSO Pipeline
This pipeline was commissioned to take input from the Level 2 pipeline -- having been processed through the NIRCam `ncdhas` pipeline -- and now being further processed through from a stack of images into a time series

**NEW METHOD**

1. Develop a single routine that inputs 
    1. String (fits file name) or array (loaded fits file)
    2. The expected location of the star (center of frame is default)
    3. Subframe size (for better center fitting)
    4. List of aperture radii (or a float for a single aperture radii)
2. This routine will load a single fits file or list of fits files (one at a time; recursive?)
3. For each single, or recursively for a list of fits files, 
    1. load the data.
    2. Computer the time element
    3. subtract the background (store background level)
    4. isolate the star into a subframe
    5. Cross-correlate a Gaussian (or JWST psf) with the image to find predicted center (store CC center)
    6. Gaussian fit to subframe, starting at CC center (store GS center, width, amplitude)
    7. Perform apeture photometry with each radius given at the beginning (store aperture radii as a function of radius)

This routine ensures that the user can manipulate the inputs as needed. Users can either send a single fits array, a set of fits array, a single string with the location of a fits file, or a list of strings with the location of several fits files.

The result will be a 'DataFrame' of the same depth as the input structure, containing (labeled as keys) the 
- 'sky background'
- 'cross correlation center'
- 'gaussian center'
- 'gaussian width'
- 'gaussian ampitude'
- 'aperture photometry dictionary' or 'aperture photometry dataframe'
    - the keys to the aperture photometry dictionary or data frame will be the float values of the aperture radii
- 'time' (in days?)
