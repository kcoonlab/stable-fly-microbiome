## Overview 
**Raw data and scripts for:**
Stable flies are bona fide carriers of mastitis-associated bacteria

## Authors 
* Andrew J. Sommer
* Julia E. Kettner
* Kerri L. Coon - kerri.coon@wisc.edu

## Analysis overview 
Scripts for each analysis are written in R. Each directory contains necessary files and code to recreate each figure, table, and associated statistical analyses reported in the manuscript. To repeat the analysis, clone the repository, and then run each script. Do not `cd` into the cloned repository. 

**Example**
* Create a new project in RStudio. To run the script to recreate Fig. 1 in the manuscript: 
	* Navigate to the terminal window and `git clone https://github.com/kcoonlab/stable-fly-microbiome`.
	* Open the script `Fig1.R` from the files panel window.
	* Install required packages. 
	* `cmd enter` from line `1`.

**Before getting started**
* Run the script `phyloseq-object.R` (code from line `1` to `42`) to generate a phyloseq object from the appropriate qiime artifacts. This object will be referenced in some downstream analyses.

**Acknowledgements**
* Portions of the following scripts are based on code produced by Jaimie R. West and Thea Whitman at the University of Wisconsin-Madison (see https://github.com/jaimiewest/Soil-Mixing and https://doi.org/10.1093/femsec/fiac112 for more information).

## Recreate the manuscript figures and associated statistical analyses
Once the repository has been cloned (above), recreate each figure/table/data analysis as follows: 

**Fig. 2: Community-level ASV richness in Arlington- and DCC-derived samples, by sample source**
* Left panel (Arlington samples) - Script: `Fig2.R`: run code from line `1` to `141`
* Right panel (DCC samples) - Script: `Fig2.R`: run code from line `1` to `154`

**Fig. 3: Relative abundance of bacterial orders in Arlington- and DCC-derived fly and manure samples, by sampling date**
* Script: `Fig3.R`: run code from line `1` to `90`

**Fig. 4: PCoA of Bray-Curtis dissimilarities of community relative abundances, colored by sample source**
* Left panel (Arlington samples) - Script: `Fig4.R`: run code from line `1` to `75`
* Right panel (DCC samples) - Script: `Fig4.R`: run code from line `1` to `95`

**Fig. 5: Bacterial ASV analysis**
* Fig 5A - Script: `Fig5.R`: run code from line `1` to `240` (Arlington samples) and from `1` to `724` (DCC samples)
* Fig 5B - Script: `Fig5.R`: run code from line `1` to `278` (Arlington samples) and from `1` to `762` (DCC samples)
* Fig 5C - Script: `Fig5.R`: run code from line `1` to `498` (Arlington samples) and from `1` to `982` (DCC samples)

**Fig. 7: Relative abundance of commonly shared ASVs in Arlington- and DCC-derived samples, by sample source**
* Left panel (Arlington samples) - Script: `Fig5.R`: run code from line `1` to `111`
* Right panel (DCC samples) - Script: `Fig5.R`: run code from line `1` to `193`

#### Miscellaneous data analyses
* Script: `misc-analyses.R`: run code from line `1` to `XX`

#### Supplementary information 
* Fig. S1 - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S2A - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S2B - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S2C - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S3 - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S4A - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S4B - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S5 - Script: `supp-info.R`: run code from line `1` to `XX`
* Fig. S6 - Script: `supp-info.R`: run code from line `1` to `XX`

## Citation 
Sommer AJ, Kettner JE, Coon KL. (2024). Microbiota composition associates with mosquito productivity outcomes in belowground larval habitats. mSphere 9(7):e0033624. https://doi.org/10.1128/msphere.00336-24.

