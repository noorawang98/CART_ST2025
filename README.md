## **CART_ST2025**  
**Integrated bioinformatics analysis of spatial transcriptome and single-cell transcriptome in CAR-T treated mouse tumors**  

### **1. Unsupervised Clustering and Cell Type Annotation**  
• **Single-cell transcriptome analysis**:  
  Pipline Script: [`sc_RNA_Reanalysis_scANCI.R`](sc_RNA_Reanalysis_scANCI.R) performs clustering and cell type annotation.  
• **Spatial transcriptome (low-resolution)**:  
  Script: [`sc_st_workflow.R`](sc_st_workflow.R) performs Unsupervised clustering and cell annotation.  
• **Spatial object creation**:  
  Tool: [`CreateSpatialObj.R`](functions/CreateSpatialObj.R) manually generates spatial objects.  
• **High-resolution deconvolution**:  
  Scripts in directory: [`Deconvolution/`](Deconvolution/)  
• **H&E alignment visualization**:  
  Script: [`hefigure.R`](Plot/hefigure.R) generates spatial plots aligned with H&E-stained images.  

### **2. Deconvolution Models**  
Tools for multi-resolution spatial deconvolution using annotated cell types:  
1. **CARD**:  
   • Repository: [YMa-lab/CARD](https://github.com/YMa-lab/CARD)  
2. **DestVI**:  
   • Documentation: [DestVI Tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/spatial/DestVI_tutorial.html)  
3. **Label Transfer with scVI**:  
   • Tutorial: [scVI Basics](https://docs.scvi-tools.org/en/stable/tutorials/index.html)  

### **3. Spatial Cell Density Statistics (kNN-based)**  
• **R implementation**:  
  Script: [`functions.R`](functions/functions.R)  
  • Calculates kNN-based cell density.  
• **Python demo**:  
  Notebook: [`Spatial_kNN_distance.ipynb`](functions/Spatial_kNN_distance.ipynb)  
  • Example workflow for kNN distance statistics.  
• **Tumor boundary identification**:  
  Script: [`Bg_ST_boundries.R`](Scoring/Bg_ST_boundries.R)  

### **4. Functional Analysis and Clinical Validation**  
#### **Trajectory Analysis and Scoring**  
• **Script directory**: [`Scoring/`](Scoring/)  
  • Cell trajectory analysis.  
  • GSVA score reclustering for CAR-T-infiltrated pixels.  

#### **TCGA Cohort Validation**  
• **CIBERSORTx Pipeline**:  
  1. **Signature Matrix Generation**:  
     Web tool: [CIBERSORTx Server](https://cibersortx.stanford.edu/)  
  2. **Human-to-Mouse Gene ID Conversion**:  
     Script: [`Bg_tcga_gene2id.R`](Clinical/Bg_tcga_gene2id.R)  
  3. **Input Preparation for CIBERSORT**:  
     Script: [`Bg_prepare_input_for_cibersortx.R`](Clinical/Bg_prepare_input_for_cibersortx.R)  
  4. **TCGA Cohort Analysis**:  
     Script: [`Bg_CIBERSORT_tcga_st_nr.R`](Clinical/Bg_CIBERSORT_tcga_st_nr.R)  
     ◦ Runs CIBERSORT across TCGA cohorts.  
     ◦ Plots survival results by cancer type.  
• **Data Sources**:  
  • Expression/meta files: Downloaded from [TIMER2.0](https://cistrome.shinyapps.io/timer/).  
