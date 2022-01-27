# Particle Size Distribution Function (PSD)
Some knowledges on cloud physics in remote sensing

# 0 Why gamma function

There are multiple ways to describe the psd in RTM. 

Bin mode

Log normal mode

Gamma mode

"The aerosol optical properties can be directly calculated based on a radiative transfer model for each bin of the CARMA size distribution. However, the Mie calculation in the current LP aerosol code requires an analytic aerosol mode, rather than bin data. In this study, we use an analyt- ical model of aerosol particle size distribution which deals with the ASD as a mean of size spectrum to accurately fit a cumulative distribution function (CDF) on the binned data using Deshler’s method (Deshler et al., 2003)." from Chen et al., 2018

Why it is called the gamma function can be found at Page 203 from Petty's Book.


# 1 The psd in OMPS LP for level 2 product:

This psd is designed for aerosol.

## 1.1 Assumptions:
1. spherical Mie scattering particles of refractive index 1.448 + 0i to calculate the effect of aerosols on LP radiance.

2. The assumed aerosol size distribution (ASD) does not vary with time or location

3. The function is described by the gamma function formulation.

## 1.2 Equation & Properties:

![Screen Shot 2022-01-26 at 3 38 37 PM](https://user-images.githubusercontent.com/52504365/151251296-3f596bca-b41a-4ace-90c4-ca5668547153.png)

n(r) is the number of particles N (r) per unit volume with a size between radius r and r+dr (cm-3 * μm-1), N0 is the total number density of aerosols (cm-3), α and β (μm-1) are the fitting parameters, and Г is Euler’s Gamma function.

The gamma function has only two parameters to be specified, the shape parameter α and the scale parameter β. α = 1.8, β = 20.5

The PSD has a relationship to effective radius

![Screen Shot 2022-01-26 at 4 01 22 PM](https://user-images.githubusercontent.com/52504365/151254305-a02c9fc1-87a5-42c5-a713-b9076ee5c010.png)


## 1.3 Measurement
Sampled for specific conditions (June-July-August average, 41°N, 20 km) corresponding to long-term balloon particle size measurements.

## 1.4 History
This psd is for LP V1.5 (I think V2.0 keeps the same).

The old LP V1.0 psd is a bimodal lognormal size distribution (BD).

The details and limitations for the old psd function can be found at Section 3 at Chen et al. (2018) 

The current V1.5/V2.0 is from Chen et al. (2018)

![fig_lpv2 0](https://user-images.githubusercontent.com/52504365/151255247-948f7e3e-18be-4581-91c4-349f04fba94a.png)


## 1.5 Code

```python
import math
import numpy as np

alpha = 1.8
beta = 20.5
num = 100

r = np.logspace(-2., 0., num=num)
n = np.zeros(num)

for i in range(num):
    tmp = (beta ** alpha) * (r[i] ** (alpha - 1))
    n[i] = tmp / math.gamma(alpha) * math.exp(-1 * r[i] * beta)
```

# 2 The psd in MODIS C6 for level 2 product:

This psd is designed for both liquid and ice clouds.

## 2.1 Equationv & Properties

The single scattering properties of liquid water clouds are calculated from Mie theory, and are integrated over a Modified Gamma droplet size distribution,

![Screen Shot 2022-01-26 at 4 15 36 PM](https://user-images.githubusercontent.com/52504365/151256131-8b6c82be-d7cb-4b09-9104-6d3b93616a5f.png)

re: Cloud Effective particle Radius

𝓋e: effective variance 𝓋e=0.10

Details see equation (2.9.1) at MODIS 06 ATBD.

So, there are

alpha - 1 = (1 - 3 * ve) / ve
beta = re * ve

The alpha = 8.

For MODIS, Re = Reff. However, this is only for spherical ice particles. Details see Page 135 at Wendish and Yang Book.

![fig_mod06](https://user-images.githubusercontent.com/52504365/151266932-5350762b-eba7-47a4-8b8f-36250c0d81a7.png)

This figure assmes re = 30 um.

## 2.2 Code

```python
import math
import numpy as np

r = np.logspace(0., 2., num=num)
n = np.zeros(num)

ve = 0.1
re = 30
num = 100

for i in range(num):
    tmp = r[i] ** ((1 - 3 * ve) / ve)
    n[i] = tmp * math.exp(-1 * r[i] / re / ve)

```

# References
Suomi-NPP OMPS LP L2 AER Daily Product ATBD, Version 2.0, Last Revised 07/15/2020

MODIS Cloud Optical Properties: User Guide for the Collection 6 Level-2 MOD06/MYD06 Product and Associated Level-3 Datasets

Chen et al. 2018, Improvement of stratospheric aerosol extinction retrieval from OMPS/LP using a new aerosol model, AMT

Petty's Book

Wendish and Yang Book

