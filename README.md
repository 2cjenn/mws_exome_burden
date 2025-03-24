# Analysis Scripts

## Develop the breast cancer [phenotype](https://github.com/2cjenn/mws_exome_burden/tree/main/Scripts/Phenotypes)

- **[LoadData.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/Phenotypes/LoadData.Rmd)** - Read in the phenotypic data and generate breast cancer phenotype, as well as other phenotypes for supplementary analyses

## [Annotation](https://github.com/2cjenn/mws_exome_burden/tree/main/Scripts/Annotation) of the WES data

1. **[Annotation.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/Annotation/Annotation.Rmd)** - Runs VEP with some plugins to annotate the data
2. **[Helix_MWS.R](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/Annotation/Helix_MWS.R)** - Scrapes annotations from the Helix API
3. **[BCAC_Annotation.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/Annotation/BCAC_Annotation.Rmd)** - Uses the annotations to create masks, following the approach of [Wilcox et al](https://www.nature.com/articles/s41588-023-01466-z)

## [Burden](https://github.com/2cjenn/mws_exome_burden/tree/main/Scripts/Burden) tests

- **[RunRegenie.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/Burden/RunRegenie.Rmd)** - Launches Regenie association tests

## Conduct the [meta-analysis](https://github.com/2cjenn/mws_exome_burden/tree/main/Scripts/MetaAnalysis)

- Read in results, conduct meta-analysis and produce output tables
  -  **[genebased_meta.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/MetaAnalysis/genebased_meta.Rmd)**
- Code for additional analyses
  - Meta-analysis including FinnGen
    - **[FinnGen.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/MetaAnalysis/FinnGen.Rmd)**
  - Family history by proxy, as per [Wilcox et al](https://www.nature.com/articles/s41588-023-01466-z)
    - **[FaH_Proxy.R](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/MetaAnalysis/FaH_Proxy.R)** and **[FaH_Proxy.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/MetaAnalysis/FaH_Proxy.Rmd)**
   
## [Plot](https://github.com/2cjenn/mws_exome_burden/tree/main/Scripts/Plotting) Supplementary Figure 1

- **[PanelPlot.Rmd](https://github.com/2cjenn/mws_exome_burden/blob/main/Scripts/Plotting/PanelPlot.Rmd)** - Create the panel figure combining forest plots for genes significant in the meta-analysis
