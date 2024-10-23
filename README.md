# Stochastic skeleton model for the Madden-Julian Oscillation with time-dependent observation-based forcing

This repository implements and analyzes the Stochastic MJO Skeleton Model with time-dependent observation-based forcing functions. This work corresponds to the results presented in Chapter 9 of my thesis: [Complex Systems Perspectives on Large-Scale Weather and Climate Variability Patterns, Noémie Ehstand](https://www.researchgate.net/publication/382743987_Complex_systems_perspectives_on_large-scale_weather_and_climate_variability_patterns).

## Running the notebooks

To run the notebooks in this repository, ensure you have Python and Julia installed along with the necessary packages (listed below). You can clone the repository using:

```bash
git clone https://github.com/noemiee/skeleton-model.git
cd skeleton-model
```


Required Julia packages: `LinearAlgebra`, `FFTW`, `PyPlot`, `Statistics`, `Interpolations`, `Random`, `Distributions`

Required python packages: `numpy`, `matplotlib`, `pandas`, `pickle`, `netCDF4`, `cftime`, `datetime`, `xarray`, `scipy`

## Data Sources
- NCEP-NCAR reanalysis latent heat net flux (Kalnay et al., 1996) is available [here](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)
- NCEP global precipitation climatology project data (Adler et al., 2017; Huffman et al., 2001) is available [here](https://www.ncei.noaa.gov/products/climate-data-records/precipitation-gpcp-monthly).

## Notebooks Overview

**`forcing-profiles.ipynb`**

This notebook computes the model’s forcing based on the methodology presented in [Ogrosky and Stechman (2015)](https://doi.org/10.1002/qj.2552), while retaining the time-dependence of the profiles.
 
**`Stochastic-skeleton-model-main.ipynb`**

This notebook contains the implementation of the stochastic MJO skeleton model, along with post-processing and analysis of the model's outputs.

**`Stochastic-skeleton-model-plots.ipynb`**

This notebook is designed specifically to reproduce the plots presented in my thesis. The data required for the plots is available [here](https://cloud.ifisc.uib-csic.es/nextcloud/index.php/s/9j9QBPZ4ZgkkNHD).

