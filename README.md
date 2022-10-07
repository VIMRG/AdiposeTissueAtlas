# White Adipose Tissue Compositional and Transcriptional Patterns with Progressive Glucose Intolerance in Persons with HIV

Samuel S. Bailin<sup>1</sup>, Jonathan A. Kropski<sup>2,3,4</sup>, Rama D. Gangula<sup>5</sup>, LaToya Hannah<sup>6</sup>, Joshua D. Simmons<sup>6</sup>, Mona Mashayekhi<sup>6</sup>, Fei Ye<sup>7</sup>, Run Fan<sup>7</sup>, Abha Chopra<sup>8</sup>, Ramesh Ram<sup>8</sup>, Simon Mallal<sup>1,5,8-10</sup>, Christopher M. Warren<sup>5</sup>, Spyros Kalams<sup>1,5</sup>, Curtis L. Gabriel<sup>11</sup>, Celestine N. Wanjalla<sup>1,10</sup>, John R. Koethe<sup>1,3,10</sup>

<sup>1</sup>Department of Medicine, Division of Infectious Diseases, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>2</sup>Department of Medicine, Division of Allergy and Pulmonology, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>3</sup>Veterans Affairs Tennessee Valley Healthcare System, Nashville, Tennessee, USA
<sup>4</sup>Deparment of Cell and Developmental Biology, Vanderbilt University, Nashville, Tennessee, USA
<sup>5</sup>Tennessee Center for AIDS Research, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>6</sup>Department of Medicine, Division of Diabetes, Endocrinology, and Metabolism, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>7</sup>Department of Biostatics, Vanderbilt University, Nashville, Tennessee, USA
<sup>8</sup>Insitute for Immunology and Infectious Diseases, Murdoch University, Perth, WA, Australia
<sup>9</sup>Vanderbilt Technologies for Advanced Genomics, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>10</sup>Center for Translational Immunology and Infectious Diseases, Vanderbilt University Medical Center, Nashville, Tennessee, USA
<sup>11</sup>Department of Medicine, Division of Gastroenterology, Hepatology, and Nutrition, Nashville, Tennessee, USA


## **Abstract**
Subcutaneous adipose tissue (SAT) is a critical regulator of systemic metabolic homeostasis. Persons with HIV (PWH) have an increased risk of metabolic diseases and have significant alterations in SAT immune cells compared with the general population. We generated a comprehensive SAT atlas from 59 PWH with a spectrum of metabolic diseases to characterize cellular compositional and transcriptional changes associated with glucose intolerance and demographic characteristics. In PWH, glucose intolerance was associated with increased lipid-associated macrophages, CD4+ and CD8+ T effector memory cells, and decreased perivascular macrophages. We observed a coordinated intercellular regulatory program enriched for genes related to inflammation and lipid-processing with glucose intolerance that was also present in diabetic HIV-negative persons. Finally, CD4+ tissue resident cells were most strongly associated with altered whole tissue expression of genes related to lipid and glucose metabolism. These data provide a comprehensive assessment of SAT immune cells that may contribute to glucose intolerance.

## **Description**
This repository contains the code necessary to reproduce the figures in the manuscript. Some steps are non-deterministic and therefore, may not yield the same results.

## **Links**

### **External**
[NIH GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198809)<br/>
[Interactive Website](http://vimrg.app.vumc.org/)<br/>
[Manuscript](https://somewebsite.com)<br/>

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
| [Pseudotime](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/Pseudotime_PVM.R) |
| [Single Cell DGE](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/SingleCellDGE.R) |
| [Dialogue](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/Analyses/Dialogue.R) |

### **HIV-negative Integration Code**
| Scripts |
| --- |
| [Initial Integration](https://github.com/VIMRG/AdiposeTissueAtlas/blob/main/Scripts/HIVnegative Analyses/Initial_Merge_HIVnegative.R) |

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



