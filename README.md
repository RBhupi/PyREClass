# Note
This package is merged into DoE ARM Radar Toolkit (PyART) package and can be found [here](https://arm-doe.github.io/pyart/API/generated/pyart.retrieve.conv_strat_raut.html).

# PyREClass

Classification of radar reflectivity data using Wavelet transform.

## Description

A new method for separating convective, intermediate, and stratiform rainfall based on the wavelet scale analysis. The method uses *à trous* wavelet transform (WT) to separate heterogeneous convective features from the relatively smooth stratiform field. The heterogeneous region is further split into a convective class for the large wavelet coefficients and an intermediate class for the small wavelet coefficients. There is an improvement in the estimation of convective-stratiform fraction and frequency distribution of rain rates in comparison to the well established Steiner method. This new method correctly reclassified a significant fraction of Steiner’s convective regions.

## Reference

Raut, B. A., Louf, V., Gayatri, K., Murugavel, P., Konwar, M., & Prabhakaran, T. (2020). A Multiresolution Technique for the Classification of Precipitation Echoes in Radar Data. *IEEE Transactions on Geoscience and Remote Sensing*, 1–7. [https://doi.org/10.1109/TGRS.2020.2965649](https://doi.org/10.1109/TGRS.2020.2965649)
