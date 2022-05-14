# Directories

## Main Analysis directories
- `in-vitro-pipeline/`: analysis of the smart-seq2, in-vitro cytokine-stimulation data
- `10xAll/`: analysis of the in-vivo, 10x, single-cell data
- `cytokine_rnaseq/`: analysis of the bulk rna-seq from in-vitro cytokine stimulation experiments
- `cd160_knockout/`: analysis of the ex-vivo T cell rna-seq from cd160-knockout experiment

See README files inside each subdirectory above for further instructions

## Other utility directories

- `GWAS/`: downloaded GWAS associations from delange/huang publications
- `common/`: python package defining consistent names/colors for figures
- `packages/`: additional custom python packages used for analysis/visualization
- `pipelineScripts/`: python/R scripts used in Snakemake pipelines
- `Signatures/`: signatures used in various enrichment analyses
- `SupplementaryFiles/`: Generation of supplementary tables for the manuscript

# Outputs
## Figure Panels with custom code

- Figure 1E - tSNE plots for in-vitro single-cell cytokine experiments
    - /in-vitro-pipeline/Th1_Pipeline/scripts/MakeFiguresTsne.py
- Figure 1F - Violin plots - GFP associated differences in-vitro single-cell
    - /in-vitro-pipeline/Th1_Pipeline/scripts/plotViolins.py
- Figure 1G - Extent of GFP-associated differential expression, WT vs. IL23R KO
    - /in-vitro-pipeline/Th1_Pipeline/scripts/MakeEffectSize.py
- Figure 1H - Genotype-dependent GFP-associated differential expression
    - /in-vitro-pipeline/Th1_Pipeline/scripts/Interaction_Volcano.py
<br>

- Figure 2A - PCA of in-vitro cytokine experiment - bulk rna-seq
    - /cytokine_rnaseq/Scripts/plot_pca.py
- Figure 2B - Heatmap of in-vitro cytokine experiment - bulk rna-seq
    - /cytokine_rnaseq/Scripts/plot_heatmap.py
- Figure 2C - Heatmaps of single-cell in-vitro cytokine rna-seq
    - /in-vitro-pipeline/Figures/plot_heatmaps.py
<br>

- Figure 4C - In-vivo single-cell UMAPs of ex-vivo T cells
    - /10xAll/Figures/combined_all/tissue_umap.py
- Figure 4D - Tissue-associated gene expression
    - /10xAll/Figures/combined_all/tissue_de.py
- Figure 4E - In-vivo single-cell UMAPs
    - /10xAll/Figures/combined_all/genotype_umap.py
- Figure 4F - Genotype-associated differential expression by tissue
    - /10xAll/Figures/combined_all/tissue_geno.py
<br>

- Figure 5A/B - Colon LPL UMAPs, Genotype and Cluster
    - /10xAll/Figures/LPL_combined/plotUMAP.py
- Figure 5C - Cluster proportions and prevalence by Genotype
    - /10xAll/Figures/LPL_combined/plotClusterProportions.py
- Figure 5D - Top cluster associated genes (LPL-2, LPL-9, and LPL-7)
    - /10xAll/Figures/LPL_combined/plotClusterDot.py
- Figure 5E - Cluster DE of IBD GWAS genes
    - /10xAll/Figures/LPL_combined/plotGWASHeatmap.py
<br>

- Figure 6A - Prioritized ranking of IL-23R dependent drivers of Th1-mediated colitis
    - /10xAll/Figures/Targets/makeFigure.py
<br>

- Figure 7A - Differential genes CD160-/- vs Ctrl
    - /cd160_knockout/Scripts/make_figures.py
- Figure 7B - GSEA enrichments
    - /cd160_knockout/Scripts/make_figures.py
- Figure 7E - Differential genes heatmap
    - /cd160_knockout/Scripts/make_figures.py
- Figure 7F - single-cell cluster signatures scored on CD160-/- rna-seq
    - /cd160_knockout/Scripts/make_signature_figures.py
- Figure 7G - CD160-/- signature scored on single-cell clusters
    - /cd160_knockout/Scripts/make_signature_figures.py
<br>

- Supplementary Figure 2 - In-vitro expression of T cell subtype markers by sample
    - /in-vitro-pipeline/Th1_Pipeline/scripts/MakeFiguresTsne.py
<br>

- Supplementary Figure 5A - T cell and NK markers within colon LPL, UMAPs
    - /10xAll/Figures/LPL_combined/plotUMAPGrids.py
- Supplementary Figure 5B - In-vivo colon LPL UMAPs - genes associated with LPL-8
    - /10xAll/Figures/LPL_combined/plotUMAPGrids.py
<br>

- Supplementary Figure 6A - In-vivo colon LPL UMAPs - genes associated with LPL-7
    - /10xAll/Figures/Tr1/plotUMAPGrid.py
- Supplementary Figure 6B - Dot plot - genes associated with LPL-7
    - /10xAll/Figures/Tr1/marker_dots.py
- Supplementary Figure 6C - Tr1 signature score, LPL-7 vs remainder
    - /10xAll/Figures/Tr1/plotTr1Sig.py
- Supplementary Figure 6D - Genotype-dependent expression within LPL-7
    - /10xAll/Figures/Tr1/de_genotype.py
<br>

## Supplementary Files

- Supplementary Table 1 - In-vitro Smart-Seq2 cell meta-data and differential expression
    - SupplementaryFiles/InVitroSupplement.xlsx
- Supplementary Table 2 - In-vivo 10x cell meta-data and differential expression (all cells, between-tissue DE)
    - SupplementaryFiles/InVivoAllSupplement.xlsx
- Supplementary Table 3 - In-vivo 10x cell meta-data and differential expression (Colon and Small Intestine LPL cells, between-cluster DE)
    - SupplementaryFiles/InVivoLPLSupplement.xlsx
- Supplementary Table 4 - Details on ranking procedure for prioritized follow up of targets
    - SupplementaryFiles/RankingSupplement.xlsx
- Supplementary Table 5 - Differential expression results from bulk rna-seq on in-vitro cytokine-stimulated CD4 T cells
    - SupplementaryFiles/CytokineDESupplement.xlsx
- Supplementary Table 6 - Differential expression results from ex-vivo colonic CD160-/- vs WT CD4 T cells
    - SupplementaryFiles/CD160KOSupplement.xlsx

# Prerequisites

In order to run the analysis, the after downloading, the following needs to be set up:

1. The python package in /common must by installed
    - navigate to /common and run `pip install -e .`
2. The python packages in /packages/bio_utils and /packages/gene_enrich must be installed
    - navigate to /packages/bio_utils and run `pip install -e .`
    - navigate to /packages/gene_enrich and run `pip install -e .`

# Environment

This code was run using Python v3.6.8, R version 3.5.1
For python, conda environment specification found in `conda.yml`
For R, list of installed R packages/versions found in renv.txt
