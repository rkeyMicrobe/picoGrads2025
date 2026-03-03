

# picoGrads2025

**Author:** Rebecca Key (R.S.Key; ORCiD 0000-0002-9516-1645)  
**Repository:** [picoGrads2025](https://github.com/rkeyMicrobe/picoGrads2025)  
**Manuscript:** Submitted to *mSystems* on May 29, 2025. Accepted December 21, 2025. Published January 16, 2026. 

---

This repository contains the analysis code and supporting files for:  
**_Picophytoplankton Implicated in Productivity and Biogeochemistry in the North Pacific Transition Zone_**
**Publication:** [mSystems DOI](https://doi.org/10.1128/msystems.00801-25)
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

This repository contains 12 R scripts. It is a step-by-step pipeline (1 to 12) for the entire analysis. Basic amplicon analysis using phyloseq, multivariate mixed linear modeling using sommer, and network analyses using WGCNA and SpiecEasi.

**How to run:**

- Start with `01_qiime2.R`

- Proceed step-by-step through each numbered script, until the last step: `12_speicEasiAnalysis.R`

⚠️ For`10_wgcnaPowerTest.R` — re-run this script for each power value you wish to test (power = ??). I didn't make it a loop script so you have to repetively run the script after changing the power variable.

Simple and reproducible: As long as you have the starting inputs for `01_qiime2.R`, you can generate all dataframes, figures, and tables used in the manuscript.

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
- Processed dataframes: [Zenodo Record](https://zenodo.org/records/17288024)  
- Raw amplicon dataframes: [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1302492/)
- Or contact the authors for `01_qiime2.R` input files.

This repository contains scripts only. **The total size after running everything is 5.1GB!** 

## Acknowledgments

We thank the scientific team and crew of the R/V Kaʻimikai-O-Kanaloa (KOK1606; Gradients 1), R/V Marcus G. Langseth (MGL1704; Gradients 2), and R/V Kilo Moana (KM1906; Gradients 3) and the operational staff of the Simons Collaboration on Ocean Processes and Ecology (SCOPE) team. We also thank Bennet Lambert for early assistance with ASV data processing.
This work was supported by grants from the Simons Foundation (Awards 823165 and 999397 to B.P.D.; Award 721244 to E.V.A.; Award 724220 to J.P.Z.; Award 426570SP to E.V.A. and J.P.Z.; Award 00012203 to S.N.C.) as part of the SCOPE Program.

R.S.K., B.P.D., and E.V.A. conceived the research project. E.V.A., J.P.Z., M.R.G., H.F., and B.P.D. planned fieldwork sampling design. R.S.K. led data analysis, data interpretation, figure generation, and manuscript writing, with supervision from B.P.D, and S.N.C. and E.V.A. assisted with refinement of research design. M.R.G. collected seawater samples. R.L.M., H.F., and M.R.G. prepared samples for ASV sequencing. All authors contributed to reviewing and revising the manuscript.

## Contact

For questions, please reach out to **Rebecca S. Key** (rkeyMicrobe [at] proton [dot] me), **Bryndan P. Durham** (b [dot] durham [at] ufl [dot] edu) or open a [GitHub Issue](https://github.com/rkeyMicrobe/picoGrads2025/issues).



