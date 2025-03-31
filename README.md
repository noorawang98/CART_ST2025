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
  • Cell-cell communication analysis.  
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
----
### **Reference**  
1. **Key Methodology Citations**:  
   • **CARD**:  
     ```plaintext
     Ma, Y., & Zhou, X. (2021). Spatially informed cell-type deconvolution for spatial transcriptomics. 
     Nature Biotechnology. https://doi.org/10.1038/s41587-021-01070-8

     ```  
   • **DestVI/scVI**:  
     ```plaintext
     Gayoso, A., et al. (2022). Joint probabilistic modeling of single-cell and spatial transcriptomes 
     with scvi-tools. Nature Methods. https://doi.org/10.1038/s41592-021-01326-x

     ```  
   • **CIBERSORTx**:  
     ```plaintext
     Newman, A.M., et al. (2019). Determining cell type abundance and expression from bulk tissues with 
     digital cytometry. Nature Biotechnology. https://doi.org/10.1038/s41587-019-0114-2

     ```  
   • **Single-cell analysis (Seurat/Scanpy)**:  
     ```plaintext
     Satija, R., et al. (2015). Spatial reconstruction of single-cell gene expression data. 
     Nature Biotechnology. https://doi.org/10.1038/nbt.3192

     Wolf, F.A., et al. (2018). SCANPY: Large-scale single-cell gene expression data analysis. 
     Genome Biology. https://doi.org/10.1186/s13059-017-1382-0

     ```  

1. **Data Sources**:  
   • TCGA expression/meta files:  
     ```plaintext
     Li, T., et al. (2020). TIMER2.0 for analysis of tumor-infiltrating immune cells. 
     Nucleic Acids Research. https://doi.org/10.1093/nar/gkaa407

     ```  

---
### **Citation**
if you use our in-house scripts,please cite DOI(https://github.com/marramWang/CART_ST2025/edit/main/README.md). and processed data is available at [Spatial multi-omic profiling identifies a stroma-immune barrier driven by CCR1⁺ myeloid impeding CAR-T therapy efficacy against tumors](https://figshare.com/account/home#/collections/7744679)
