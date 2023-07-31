# Temporal dynamics of normalization reweighting (code)

These materials are a computationally reproducible version of the paper:

Baker, D.H., Marinova, D., Aveyard, R., Hargreaves, L.J., Renton, A., Castellani, R., Hall, P.,  Harmens, M., Holroyd, G., Nicholson, B., Williams, E.L., Hobson, H.M. & Wade, A.R. (2023). Temporal dynamics of normalization reweighting, BioRxiv, https://doi.org/10.1101/2023.02.10.527994

The file manuscript.Rmd is an R markdown file that will perform all analyses and figure creation, and produce a pdf version of the manuscript.

The full repository can be downloaded (cloned), and will automatically download the required packages and data files, depending on the level of analysis requested. If any data files are missing the code will attempt to download them from the OSF repository for this project:
http://doi.org/10.17605/OSF.IO/ab3yv

The 'docker' directory contains a Dockerfile and instructions for making a local computationally reproducible version of the analysis. In addition, the Docker environment is set up to run automatically on a remote server via Github Actions, each time a change is made (i.e. on a 'commit' to the repo). The output document is then posted back to the main repository (manuscript.pdf). If you want to make changes to the analysis and have these build automatically, you can fork the repository into your own account.

![autobuild](https://github.com/bakerdh/normreweight/workflows/autobuild/badge.svg)
