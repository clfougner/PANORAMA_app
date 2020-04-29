# PANORAMA web application
This repository contains the source code for the PANORAMA website hosted at https://pancancer.app. For an explanation of the web application, please the the "About" tab of the website.

This repository is associated with the following article:

Fougner, C., Höglander, E.K., Lien, T.G., Sørlie, T.S., Nord, S. & Lingjærde, O.C. A pan-cancer atlas of transcriptional dependence on DNA methylation and copy number aberrations. _bioRxiv_ (2020).

All code used in the web application can be found in the file `app.R`. The file `AllDat.rds` contains model data used in the pan-cancer view. See the above article, and https://github.com/clfougner/PANORAMA for details of how these models were generated. The remaining files in the `/data/` folder are reference files used to f.ex. identify karyotype bands and gene positions. The data used in the single model view are too large to be hosted publicly, but are based on data queried from [The Cancer Genome Atlas](https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin) and processed as described in https://github.com/clfougner/PANORAMA. Further processing of data beyond that described is only for performance purposes (i.e. splitting genome-wide files into per-gene files).
