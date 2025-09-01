# Workshop on py3SEB modelling framework
Digital *Jupyter Book* collection developed for the **Workshop on py3SEB modelling framework**.

## Installation
You can use [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/VicenteBurchard/py3seb-workshop/master) for runing the notebooks in the cloud or...

You should have these programs installed:

* Python and/or [Anaconda](https://www.anaconda.com/download/success). 
* [Git](https://git-scm.com/downloads)
* [QGIS](https://qgis.org/en/site/forusers/download.html) (not necessary but recommended to visualize outputs)

We need to install all the library requisites for this tutorial included in [requisites.txt](./requirements.txt) or in [environment.yml](./environment.yml) 

Install the requirements either with conda/mamba:

`mamba env create -f environment.yml`

or

`conda env create -f environment.yml`

or with pip:

`pip install -r requirements.txt`


## Run the interactive book
Run the following command to initialize the book:

`bash run_workshop.sh`

or 

`run_workshop.bat` if you are working under Windows.

Open your web browser to [`http://localhost:3000`](http://localhost:3000)

## Contents
## Day 1 (Introduction and proximal sensing)

1. Physical principles (2 hours).
    
    a. [Radiation transfer](./101-Net_radiation.ipynb)

    b. [Turbulent exchange of heat and vapour](./102-Turbulence_and_sensible_heat_flux.ipynb)
 
2. Introduction to TSEB and 3SEB (2 hours)
    
    a. [How TSEB works?: a quick overview of low-level code](./103-TSEB_introduction.ipynb)
	
	b. [How 3SEB works?: a quick overview of low-level code](./104-3SEB_introduction.ipynb)

## Day 2 (UAV high resolution imagery)

3. Running pyTSEB and py3SEB with UAV imagery (4 hours)
    
    a. [Retrieving canopy and soil temperatures with UAV](./201-UAV_canopy_and_soil_temperatures.ipynb)

    b. [Running TSEB-2T](./202-UAV_TSEB-2T.ipynb)

    c. [Running 3SEB with high resolution imagery](./203-UAV_3SEB.ipynb)
    

        
## Day 3 (Satellite imagery)
  
4. Satellite-based implementation (4 hours)

    a. [Copernicus TSEB and TSEB](./501-Copernicus_TSEB_3SEB.ipynb)
	
	b. [Retrieval of biophysical traits with RTM](./302-Biophysical_Traits_RTM.ipynb)


## Bonus Material

5. other tools, applications and developments
  
    a. [TSEB in row crops](./B01-Row_crops.ipynb)

    b. [Shuttleworth-Wallace TSEB](./B02-TSEB-SW.ipynb)
        
    c. [Automatic co-registration between mosaics](./B03-Mosaics_corregistration.ipynb)

    d. [Derivation of stomata conductance and link to CO$_2$ fluxes](./B04-stomata_conductance.ipynb)


## License
Creative Commons Attribution-ShareAlike 4.0 International.

This work is licensed under Attribution-ShareAlike 4.0 International. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/

This license requires that reusers give credit to the creator. It allows reusers to distribute, remix, adapt, and build upon the material in any medium or format, even for commercial purposes. If others remix, adapt, or build upon the material, they must license the modified material under identical terms.

  - BY: Credit must be given to you, the creator.
  - SA: Adaptations must be shared under the same terms. 
  



