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

---
Load All Necessary Libraries and Functions
---

    `pylab`      : combination of array manipulation and plotting functions
    `matplotlib` : specialized plotting functions
    `numpy`      : array more manipulation functions
    `pandas`     : dataframe -- more advanced array / table -- functions
    `photutils`  : astropy associated package for aperture photometry
    `astropy`    : `modeling` : access linear and gaussian functions with astropy formatting
                   `fitting`  : access to astropy fitting routines
    `glob`       : grab list of files in directory

    -- Not Used Yet --
    `astroML`    : better histogram function for plotting
    `sklearn`    : `externals`: imports operating system (storage) level function (i.e. joblib)
    `statsmodels`: `robust`   : robust statistical modeling packages; `scale.mad` == median average distance
    `sys`        : python-os level functions (i.e. path)
    `time`       : compute and convert current timestamps from python / os

```python
# Matplotlib for Plotting
%matplotlib inline
from pylab              import gcf, sort, linspace, indices, std, empty, concatenate, pi, sqrt, ones, diag, inf
from pylab              import rcParams, array, get_current_fig_manager, twinx, figure, subplots_adjust

from matplotlib.ticker  import MaxNLocator
from matplotlib         import style
from matplotlib         import pyplot as plt

# Numpy & Pandas for Array and DataFrame Manipulation
from numpy              import min, max, median, mean, zeros, empty
from numpy              import ones, where, arange, indices
from pandas             import DataFrame, read_csv, scatter_matrix

# Astropy for Aperture Photometry and Fits read/write
from photutils          import CircularAperture, aperture_photometry
# from astroML.plotting   import hist
from astropy.modeling   import models, fitting
from astropy.io         import fits

# Built in Libraries for directory
from glob               import glob

# Adam Ginsburg
from image_registration import cross_correlation_shifts

# Data Storage from Sci-kits
# from sklearn.externals import joblib

# Style Components
# from seaborn            import *

# from socket             import gethostname
# from statsmodels.robust import scale
# from sys                import exit, stdout
# from time               import time

style.use('fivethirtyeight')
```

This is an example input for the requests below. The directory contains JWST-NIRCam fits files within it
    - only works on Jonathan Fraine's Laptop
    - soon to 'upgrade' to working on surtr

'/Users/jonathan/Research/NIRCam/CV3/StabilityTest/fitsfilesonly/reduced_orig_flags/redfits/NRCN821CLRSUB1-6012172256_1_481_SE_2016-01-12T18h00m43.red/'

There is also a test file in the current working directory named `'fits_input_file.txt'`. It was creating using the bash 'script'

```bash
cd /Users/jonathan/Research/NIRCam/CV3/StabilityTest/fitsfilesonly/reduced_orig_flags/redfits/NRCN821CLRSUB1-6012172256_1_481_SE_2016-01-12T18h00m43.red/

ls > fits_input_file.txt
```

---
Responding to the inquiry with (including appostraphes) either 

`'fits_input_file.txt'` 

or 

`'/Users/jonathan/Research/NIRCam/CV3/StabilityTest/fitsfilesonly/reduced_orig_flags/redfits/NRCN821CLRSUB1-6012172256_1_481_SE_2016-01-12T18h00m43.red/'` 

is successful
---

Request Directory with a Set of Fits Files OR a Text File with the Same List
---
```python
list_of_data_file_types = ['.txt', '.dat', '.csv']
nircam_data = DataFrame()
found       = False
DataDir     = input()

for filetype in list_of_data_file_types:
    if filetype in DataDir:
        nircam_data['fitsfilenames'] = read_csv(DataDir)
        found = True

if not found:
    nircam_data['fitsfilenames'] = glob(DataDir+'/*')
```

Compute Julian Data from Header
---

This function is a wrapper for `julian_date` in the `jd.py` package (soon to be converted to `julian_date.py` package.
It's utility is in taking in the time stamps from the headers and converting them to the julian date; to be saved in the 'master' data frame below.

```python
def get_julian_date_from_gregorian_date(*date):
    """gd2jd.py converts a UT Gregorian date to Julian date.
    
    Functions for JD <-> GD conversion, 
      courtesy of Ian Crossfield at 
      http://www.astro.ucla.edu/~ianc/python/_modules/date.html
    
    Downloaded from Marshall Perrin Github at
        https://github.com/mperrin/misc_astro/blob/master/idlastro_ports/gd2jd.py
    
    Usage: gd2jd.py (2009, 02, 25, 01, 59, 59)

    To get the current Julian date:
        import time
        gd2jd(time.gmtime())

    Hours, minutes and/or seconds can be omitted -- if so, they are
    assumed to be zero.

    Year and month are converted to type INT, but all others can be
    type FLOAT (standard practice would suggest only the final element
    of the date should be float)
    """
    verbose=False
    if verbose: print date
    #print date[0]
    #date = date[0]

    date = list(date)
    
    if len(date)<3:
        print "You must enter a date of the form (2009, 02, 25)!"
        return -1
    elif len(date)==3:
        for ii in range(3): date.append(0)
    elif len(date)==4:
        for ii in range(2): date.append(0)
    elif len(date)==5:
        date.append(0)

    yyyy = int(date[0])
    mm = int(date[1])
    dd = float(date[2])
    hh = float(date[3])
    min = float(date[4])
    sec = float(date[5])

    UT=hh+min/60+sec/3600


    total_seconds=hh*3600+min*60+sec
    fracday=total_seconds/86400

    if (100*yyyy+mm-190002.5)>0:
        sig=1
    else:
        sig=-1

    JD = 367*yyyy - int(7*(yyyy+int((mm+9)/12))/4) + int(275*mm/9) + dd + 1721013.5 + UT/24 - 0.5*sig +0.5

    months=["January", "February", "March", "April", "May", "June", "July", "August", 
                "September", "October", "November", "December"]

    # Now calculate the fractional year. Do we have a leap year?
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]
    if (yyyy%4 != 0):
        days=daylist2
    elif (yyyy%400 == 0):
        days=daylist2
    elif (yyyy%100 == 0):
        days=daylist
    else:
        days=daylist2

    daysum=0
    for y in range(mm-1):
        daysum=daysum+days[y]
    daysum=daysum+dd-1+UT/24

    if days[1]==29:
        fracyear=yyyy+daysum/366
    else:
        fracyear=yyyy+daysum/365
    if verbose: 
        print yyyy,mm,dd,hh,min,sec
        print "UT="+`UT`
        print "Fractional day: %f" % fracday
        print "\n"+months[mm-1]+" %i, %i, %i:%i:%i UT = JD %f" % (dd, yyyy, hh, min, sec, JD),
        print " = " + `fracyear`+"\n"
    # print dd,mm,yyyy, hh,min,sec, UT



    return JD
```
