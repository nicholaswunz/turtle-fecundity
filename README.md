# Sea turtle reproductive-energy output
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)


This repository contains code and data needed to reproduce the article:

**Wu N. C., Rusli M. U., Broderick A. C., & Barneche D. R.** (In review) Scaling of sea turtle reproduction reconciles theory and conservation strategies at the global scale. *Global Ecology and Biogeography*,


**When using the data or code from this project, please cite it as:**

**Wu N. C., Rusli M. U., Broderick A. C., & Barneche D. R** (2022) nicholaswunz/turtle-fecundity: Accepted version of paper data and code of manuscript: Scaling of sea turtle reproduction reconciles theory and conservation strategies at the global scale (Global Ecology and Biogeography). *Zenodo*. DOI: [![DOI](https://zenodo.org/badge/201723328.svg)](https://zenodo.org/badge/latestdoi/201723328)

**Raw data**
- `meta_raw_data.csv`   - Raw data for extracted sea turtle body size and reproductive parameters used for the analysis.
- `egg_energy_data.csv` - Raw data for extracted eggs-size and egg-energy used for the analysis.
- `phylo_cor.rds`       - Phylogenetic co-variance matrix for the analysis.
- `phylo_tree.rds`      - Phylogenetic tree to used to produce Figure 1c phylogeny.

- `world_turtle.csv`    - Raw data for the percentage protected beaches by country per species used for the analysis.
- `world_turtle_grouped.csv` - Raw data for `world_turtle.csv` in long format.

- `nesting_year.csv`    - Raw data for the Redang Island turtle survey (1993-2019) used for the analysis.

**R codes**
- `species_analysis.R` - Code to clean and analyse data, and produce figures of all sea turtle species from the systematic search.
- `spatial_analysis.R` - Code to clean and analyse data, and produce Figure 2 - Global distribution of sea turtle nesting sites with current protected areas by countries.
- `green_turtle_analysis.R` - Code to clean and analyse data, and produce figures for the Redang Island turtle survey (1993-2019). 

**Extra files**
- `xxx.PDF` - Supplementary file includes statistical outcomes, additional figures, and descriptions from the main document.

## Abstract
**Aim**: Body size of marine megafauna can influence population dynamics because larger females have disproportionally greater reproductive output. We explored how this size scaling relationship can affect predictions of population-size structure in nesting sea turtles by combining a phylogenetically controlled meta-analysis with a long-term field nesting survey.  
**Location**: Global (meta-analysis), and Malaysia (field survey).  
**Time period**: Present.  
**Major taxa studied**: Sea turtles.  
**Methods**: We extracted body size and reproductive parameters of all sea turtle species from the literature and estimated the reproductive-energy output using allometric models. We then examined the relationship between body size and the proportion of protected nesting beaches by country as an indicator of conservation efforts on body size. Long-term monitoring (1993–2019) of body size and nesting data on green turtles (*Chelonia mydas*) from Redang Island (Malaysia) was used to examine temporal changes in body size and fecundity, and to test whether the size scaling of fecundity was isometric (linear) or allometric (curvilinear).  
**Results**: We show that the total reproductive-energy output of larger nesting females was disproportionately greater in all sea turtle species. We found no strong correlations for countries with higher proportion of protected nesting sites with female size. Finally, we show that scaling-derived calculations of population-level yearly reproductive output in the green turtle population from Redang Island were more accurate when using a hyperallometric (rather than an isometric) relationship at the individual level.  
**Main conclusions**: Understanding ecosystem function and conservation effort requires accurate predictions of population trends. Our findings highlight the necessity to account for scaling effects of body size in predicting anthropogenically-mediated population shifts, as well as the need to protect large females in order to facilitate effective population replenishment.

**Keywords:** allometry, body size, marine megafauna, metabolic theory, reproductive output


## License
This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).
