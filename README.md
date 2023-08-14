# Changes in Subcutaneous White Adipose Tissue Cellular Composition and Molecular Programs Underlie Glucose Intolerance in Persons with HIV

Samuel S. Bailin<sup>1</sup>, Jonathan A. Kropski<sup>2,3,4</sup>, Rama D. Gangula<sup>5</sup>, LaToya Hannah<sup>6</sup>, Joshua D. Simmons<sup>6</sup>, Mona Mashayekhi<sup>6</sup>, Fei Ye<sup>7</sup>, Run Fan<sup>7</sup>, Abha Chopra<sup>8</sup>, Ramesh Ram<sup>8</sup>, Simon Mallal<sup>1,5,8-10</sup>, Christopher M. Warren<sup>5</sup>, Spyros Kalams<sup>1,5</sup>, Curtis L. Gabriel<sup>11</sup>, Celestine N. Wanjalla<sup>1,10</sup>, John R. Koethe<sup>1,3,10</sup>

<sup>1</sup>Department of Medicine, Division of Infectious Diseases, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>2</sup>Department of Medicine, Division of Allergy, Pulmonology, and Critical Care Medicine, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>3</sup>Veterans Affairs Tennessee Valley Healthcare System, Nashville, Tennessee, USA
<sup>4</sup>Department of Cell and Developmental Biology, Vanderbilt University, Nashville, Tennessee, USA
<sup>5</sup>Tennessee Center for AIDS Research, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>6</sup>Department of Medicine, Division of Diabetes, Endocrinology, and Metabolism, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>7</sup>Department of Biostatics, Division of Epidemiology, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>8</sup>Insitute for Immunology and Infectious Diseases, Murdoch University, Perth, WA, Australia
<sup>9</sup>Vanderbilt Technologies for Advanced Genomics, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>10</sup>Center for Translational Immunology and Infectious Diseases, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>11</sup>Department of Medicine, Division of Gastroenterology, Hepatology, and Nutrition, Nashville, Tennessee, USA


## **Abstract**
Subcutaneous adipose tissue (SAT) is a critical regulator of systemic metabolic homeostasis. Persons with HIV (PWH) have an increased risk of metabolic diseases and significant alterations in the SAT immune environment compared with the general population. We generated a comprehensive single-cell multi-omic SAT atlas to characterize cellular compositional and transcriptional changes in 59 PWH across a spectrum of metabolic health. Glucose intolerance was associated with increased lipid-associated macrophages, CD4+ and CD8+ T effector memory cells, and decreased perivascular macrophages. We observed a coordinated intercellular regulatory program which enriched for genes related to inflammation and lipid-processing across multiple cell types as glucose intolerance increased. Increased CD4+ effector memory tissue-resident cells most strongly associated with altered expression of adipocyte genes critical for lipid metabolism and cellular regulation. Intercellular communication analysis demonstrated enhanced pro-inflammatory and pro-fibrotic signaling between immune cells and stromal cells in PWH with glucose intolerance compared with non-diabetic PWH. Lastly, while cell type-specific gene expression among PWH with diabetes was globally similar to HIV-negative individuals with diabetes, we observed substantially divergent intercellular communication pathways. These findings suggest a central role of tissue-resident immune cells in regulating SAT inflammation among PWH with metabolic disease, and underscore unique mechanisms that may converge to promote metabolic disease. 

## **Description**
This repository contains the code necessary to reproduce the figures in the manuscript. Some steps are non-deterministic and therefore, may not yield the same results.

## **Links**

### **External**
[NIH GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198809)<br/>
[Interactive Website](http://vimrg.app.vumc.org/)<br/>
[Manuscript](https://doi.org/10.1101/2022.03.21.484794)<br/>

### **Processing Code**
| Scripts |
| --- |
| [Pre-Processing](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Processing/PreProcessing.R) |
| [Doublet Finder](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Processing/DoubletFinder.R) |
| [Initial Integration](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Processing/Initial_Merge.R) |
| [Doublet Removal](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Processing/Doublet_Removal.R) |
| [Subset Analysis](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Processing/Subset_Analysis.R) |
| [Final Integration](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Processing/Final_Integration.R) |

### **Analyses Code**
| Scripts |
| --- |
| [Pseudobulk](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/Pseudobulk.R) |
| [Pseudotime](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/Macrophage_Pseudotime.R) |
| [Single Cell DGE](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/ClusterDGE.R) |
| [Dialogue](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/Dialogue.R) |
| [CellChat](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/CellChat.R) |

### **HIV-negative Integration Code**
| Scripts |
| --- |
| [Initial Integration](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative_Analyses/Initial_Merge_HIVnegative.R) |
| [Doublet Removal](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative_Analyses/Doublet_Removal_HIVneg.R) |
| [Subset Analysis](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative_Analyses/Subset_Analysis_HIVnegative.R) |
| [Final Integration](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative_Analyses/Final_Integration_HIVneg.R) |
| [Dialogue](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative_Analyses/Dialogue_HIVneg.R) |
| [CellChat](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative_Analyses/CellChat_HIVneg.R) |
| [Pseudobulk](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative_Analyses/Pseudobulk_HIVneg.R) |

### **Manuscript Figure Code**
| Scripts |
| --- |
| [Figure 1](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure1.R) |
| [Figure 2](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure2.R) |
| [Figure 3](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure3.R) |
| [Figure 4](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure4.R) |
| [Figure 5](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure5.R) |
| [Figure 6](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure6.R) |
| [Figure 7](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure7.R) |
| [Figure 8](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure8.R) |
| [Figure 9](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure9.R) |
| [Figure 10](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Figures/Figure10.R) |
| [Helper Functions](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Helper_Functions/Utils.R) |

### **Manuscript Supplemental Figure Code**
| Scripts |
| --- |
| [Supplemental Figure 1](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure1.R) |
| [Supplemental Figure 2](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure2.R) |
| [Supplemental Figure 3](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure3.R) |
| [Supplemental Figure 4](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure4.R) |
| [Supplemental Figure 5](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure5.R) |
| [Supplemental Figure 6](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure6.R) |
| [Supplemental Figure 7](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure7.R) |
| [Supplemental Figure 8](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure8.R) |
| [Supplemental Figure 9](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure9.R) |
| [Supplemental Figure 10](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure10.R) |
| [Supplemental Figure 11](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Supplemental_Figures/SupplementalFigure11.R) |




