---
title: Workshop on 3SEB modelling framework
subject: Tutorial
subtitle: Three-day workshop that will show the implementation of TSEB and 3SEB modelling framework using different sources of inputs 
short_title: How to 3SEB
authors:
  - name: Héctor Nieto
    affiliations:
      - Insituto de Ciencias Agrarias
      - CSIC
    orcid: 0000-0003-4250-6424
    email: hector.nieto@ica.csic.es
  - name: Benjamin Mary
    affiliations:
      - Insituto de Ciencias Agrarias
      - CSIC
    orcid: 0000-0001-7199-2885
  - name: Vicente Burchard-Levine
    affiliation:
      - Insituto de Ciencias Agrarias
      - CSIC
    orcid: 0000-0003-0222-8706
  - name: Miguel Herrezuelo
    affiliation:
      - Insituto de Ciencias Agrarias
      - CSIC
  - name: Radoslaw Guzinski
    affiliations:
      - DHI
    orcid: 0000-0003-0044-6806
license: CC-BY-SA-4.0
keywords: myst, markdown, open-science, tseb
---

# Summary
This repository contains the online material needed to run the digital notebooks (**Jupyter Notebook**) developed for the course.


# Installation
:::{tip} Prerequisites - Install Python and QGIS
:class: dropdown

You should have these programs installed:
* Python and/or [Anaconda](https://www.anaconda.com/download/success). 
* [Git](https://git-scm.com/downloads)
* [QGIS](https://qgis.org/en/site/forusers/download.html) (not necessary but recommended to visualize outputs)
:::

We need to install all the library requisites for this tutorial. 

First navigate in a terminal to the `py3SEB-workshop` folder:

(navigation)=
:::::{tab-set}
::::{tab-item} Windows
:sync: win
1. Open the `Anaconda Prompt Terminal` 
2. in the terminal navigate where to the `py3SEB-workshop` folder:
```{code} bash
cd C:\Users\<user>\Downloads\py3SEB-workshop
```
::::
::::{tab-item} Linux
:sync: linux
1. Open a Terminal (e.g. by presing `CTRL+ALT+T`)
2. in the terminal navigate where to the `py3SEB-workshop` folder:
```{code} bash
cd /home/<user>/downloads/py3SEB-workshop
```
::::
:::::


Then install the requirements either with `pip` or with `conda/mamba` (recommended):
(requirements)=
:::::{tab-set}
::::{tab-item} Pip
:sync: pip
```{code} bash
pip install requirements.txt
```
::::
::::{tab-item} Conda/Mamba
:sync: conda
```{code} bash
mamba env create -f environment
conda activate py3seb-workshop
```
::::
:::::

# Run the interactive book
:::::{tab-set}
::::{tab-item} Windows
:sync: win
1. Open the `Anaconda Prompt Terminal` 
2. in the terminal navigate where to the `py3SEB-workshop` folder:
```{code} bash
cd C:\Users\<user>\Downloads\py3SEB-workshop
```
3. Run this command
```{code} bash
run_workshop.bat
```
::::
::::{tab-item} Linux
:sync: linux
1. Open a Terminal (e.g. by presing `CTRL+ALT+T`)
2. in the terminal navigate where to the `py3SEB-workshop` folder:
```{code} bash
cd /home/<user>/downloads/py3SEB-workshop
```
3. Run this command
```{code} bash
bash run_workshop.sh
```
::::
:::::

Open your web browser to [`http://localhost:3000`](http://localhost:3000)
.

# Jupyter Book contents
:::{include} ./py3SEB_workshop_curriculum.md

:::

# License
Creative Commons Attribution-ShareAlike 4.0 International.

This work is licensed under Attribution-ShareAlike 4.0 International. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/

This license requires that reusers give credit to the creator. It allows reusers to distribute, remix, adapt, and build upon the material in any medium or format, even for commercial purposes. If others remix, adapt, or build upon the material, they must license the modified material under identical terms.

  - BY: Credit must be given to you, the creator.
  - SA: Adaptations must be shared under the same terms. 
  



