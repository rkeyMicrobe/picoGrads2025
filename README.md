

# picoGrads2025

**Author:** Rebecca Key (R.S.Key; ORCiD 0000-0002-9516-1645)  
**Repository:** [picoGrads2025](https://github.com/rkeyMicrobe/picoGrads2025)  
**Manuscript:** Submitted to *mSystems* on May 29, 2025

---

This repository contains the analysis code and supporting files for:  
**_Picophytoplankton Implicated in Productivity and Biogeochemistry in the North Pacific Transition Zone_**
**Preprint:** [bioRxiv DOI](https://doi.org/10.1101/2025.05.29.656823)

For questions or to report issues, please open an [Issue](https://github.com/rkeyMicrobe/picoGrads2025/issues) or contact **Rebecca S. Key** or **Bryndan Durham** (contact info at end)

---

## Repository Structure

```plaintext
📦 picoGrads2025
 ├── QIIME2         # This contains QIIME2 scripts (BASH) for processing RAW sequences 
 ├── data_in/       # This contains all starting inputs for the analysis
     ├── cmap       # Contains NCP, POC, and PON measurements
     ├── g1         # Gradients 1, 2016 Cruise - QIIME2 outputs
     ├── g2         # Gradients 2, 2017 Cruise - QIIME2 outputs
     ├── g3         # Gradients 3, 2019 Cruise - QIIME2 outputs
     ├── meta       # Contains sample metaData
 ├── data_out/
     ├── scriptFolder     # For each script, you have a script folder containing its outputs
         ├── dataframes   # Any RDS, feather, csv, text files will be stored here
         ├── figures      # Any svg, png will be stored here
         ├── tables       # Any tables in the form of pngs will be stored here
 ├── LICENSE
 ├── README.md            # This file will generate what you are reading right now  :>
 └── .gitignore
```

---

## Analysis Scripts

This repository includes 12 R scripts forming a custom, multi-step pipeline for amplicon sequence analysis (using phyloseq), multivariate mixed linear modeling (sommer), and network inference (WGCNA, SpiecEasi).

**How to run:**

- Start with `01_qiime2.R`

- Proceed step-by-step through each numbered script up to `12_speicEasiAnalysis.R`

⚠️ For`10_wgcnaPowerTest.R` — re-run this script for each power value you wish to test (power = ??)

Simple and reproducible: As long as the required inputs for `01_qiime2.R` are in place, you can generate all dataframes, figures, and tables used in the manuscript.


There are 12 .R scripts in total and are contained in this project folder. This is custom build multi-step pipeline for amplicon sequencing analysis using *phyloseq*, multivariate mixed linear models using *sommer*, and network inference using *WGCNA* and *SpeicEasi* packages. 

If you have any questions, please feel free to reach out!

---

## ⚙Requirements & Installation

This project was developed using **2025.08.0 Build 135** "Cucumberleaf Sunflower"

⚠️ It is recommended to run these scripts in R-Studio! 

Below is the full list of R packages used in this pipeline:

**Core tidy & plotting:**  
`tidyverse`, `lubridate`, `scales`, `reshape2`, `ggpubr`, `cowplot`, `RColorBrewer`, `viridis`, `ggrepel`, `ggforce`, `ggvenn`, `patchwork`

**Data tables & summaries:**  
`feather`, `flextable`, `gt`, `skimr`

**Ecological & statistical analyses:**  
`vegan`, `phyloseq`, `WGCNA`, `SpiecEasi`, `sommer`, `compositions`, `bestNormalize`, `igraph`, `Matrix`

**Mapping:**  
`maps`, `mapproj`

## Data Availability

Due to Github file size limits, raw and processed sequencing data are hosted externally:  
- Processed dataframes: Zenodo link Coming soon! :) 
- Raw amplicon dataframes: NCBI link coming soon!
- Or contact the authors for `01_qiime2.R` input files.

This repository contains scripts only. **The total size after running everything is 5.1GB!** 

## Acknowledgments

We thank the scientific team and crew of the R/V Kaʻimikai-O-Kanaloa (KOK1606; Gradients 1), R/V
Marcus G. Langseth (MGL1704; Gradients 2), and R/V Kilo Moana (KM1906; Gradients 3) and the
operational staff of the Simons Collaboration on Ocean Processes and Ecology (SCOPE) team. We also
thank Bennet Lambert for early assistance with ASV data processing. This work was supported by
grants from the Simons Foundation (Awards 823165 and 999397 to BPD; Award 721244 to EVA;
Award 724220 to JPZ; Award 426570SP to EVA and JPZ; Award 00012203 to SNC) as part of the
SCOPE Program.

Special thanks to the [Durham Lab](https://durhamlab.org/) and collaborators for fieldwork and discussion

## Contact

For questions, please reach out to **Rebecca S. Key** (rkeyMicrobe [at] proton [dot] me), **Bryndan P. Durham** (b [dot] durham [at] ufl [dot] edu) or open a [GitHub Issue](https://github.com/rkeyMicrobe/picoGrads2025/issues).



