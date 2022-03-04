# A Single-Cell Molecular Atlas of White Adipose Tissue Demonstrates Differences in Myeloid and Lymphoid Cell Polarization in Type 2 Diabetes and HIV Infection
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

Corresponding Author: Celestine N. Wanjalla, MD, Ph.D.;celestine.wanjalla@vumc.org

## **Abstract**
Persons with HIV (PWH) are at increased risk for cardiometabolic disease compared with HIV-negative persons. Subcutaneous white adipose tissue (scWAT) is a critical regulatory of systemic metabolic homeostasis through complex interactions with immune and stromal cells. We sought to characterize the scWAT of PWH with varying glucose tolerance and determine whether scWAT composition is different between diabetic PWH and diabetic HIV-negative persons. Using paired single-cell RNA sequencing and proteogenomic profiling, we generated a detailed molecular single-cell atlas of scWAT in humans with HIV-infection and with diabetes. We found significantly higher proportion of T cells in adipose tissue of PWH compared with HIV-negative individuals and and association between CD4 and CD8 effector memory T cells, glucose intolerance, and lipid-associated macrophages in PWH only. We further showed that pro-fibrotic stromal cells are expanded with higher proportion of lipid-associated macrophages. 

## **Description**
This repository contains the code necessary to reproduce the figures in the manuscript. Some steps are non-deterministic and therefore, may not yield the same results. We also have links to download the data and interact with our online single cell atlas.

## **Links**

### **External**
[NIH GEO](https://www.ncbi.nlm.nih.gov/geo/)<br/>
[Single-Cell Atlas]( https://imrc.shinyapps.io/shinyappmulti/)<br/>
[Manuscript](https://somewebsite.com)<br/>
[Group Website](https://koethe_lab.org)<br/>

### **Processing Code**
| Scripts | Dataset |
| --- | --- |
| [Quality Control](http://VIMRG/FATLAS/QC.rmd) |     |
| [Doublet Identification](http://VIMRG/FATLAS/SingleLane.rmd) | [Datasets](http://VIMRG/FATLAS/SingleLane.rds) |
| [Merge Data](http://VIMRG/FATLAS/Merged.rmd) | [Dataset](http://VIMRG/FATLAS/Merged.rds) |
| [Subset Analysis](http://VIMRG/FATLAS/Subset_Analysis.rmd) |    |
| [Integration](http://VIMRG/FATLAS/Integration.rmd) | [Seurat Object](http://VIMRG/FATLAS/Integrated.rds) |
| [Proportion](http://VIMRG/FATLAS/Proportion.rmd) |    |

### **Manuscript Figures**
| Scripts | Data |
| --- | --- |
| [Figure Functions](https://github.com/VIMRG/FATLAS/blob/Figures/Figure_Functions.R) |   |
| [Figure 1](https://github.com/VIMRG/FATLAS/blob/Figures/Figure1.R) | [Figure 1](http://VIMRG/FATLAS/Figure1_png.com) |
| [Figure 2](https://github.com/VIMRG/FATLAS/blob/Figures/Figure2.R) | [Figure 2](http://VIMRG/FATLAS/Figure2_png.com) |
| [Figure 3](https://github.com/VIMRG/FATLAS/blob/Figures/Figure3.R) | [Figure 3](http://VIMRG/FATLAS/Figure3_png.com) |
| [Figure 4](https://github.com/VIMRG/FATLAS/blob/Figures/Figure4.R) | [Figure 4](http://VIMRG/FATLAS/Figure4_png.com) |

### **Manuscript Supplemental Figures**
| Scripts | Data |
| --- | --- |
| [Supplemental Figure 1](https://github.com/VIMRG/FATLAS/blob/Supplemental_Figures/Supplemental_Figure1.R) | [Supplemental Figure 1](http://VIMRG/FATLAS/Supplemental_Figure1_png.com) |
| [Supplemental Figure 2](https://github.com/VIMRG/FATLAS/blob/Supplemental_Figures/Supplemental_Figure2.R) | [Supplemental Figure 2](http://VIMRG/FATLAS/Supplemental_Figure2_png.com) |
| [Supplemental Figure 3](https://github.com/VIMRG/FATLAS/blob/Supplemental_Figures/Supplemental_Figure3.R) | [Supplemental Figure 3](http://VIMRG/FATLAS/Supplemental_Figure3_png.com) |
| [Supplemental Figure 4](https://github.com/VIMRG/FATLAS/blob/Supplemental_Figures/Supplemental_Figure4.R) | [Supplemental Figure 4](http://VIMRG/FATLAS/Supplemental_Figure4_png.com) |
| [Pseudotime](http://VIMRG/FATLAS/Supplemental_Pseudotime.rmd) |    |





