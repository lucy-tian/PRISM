# PRISM: Polygenic Risk Score Transferability with Annotation-Informed Modeling

This repository contains the codebase and analysis scripts for **PRISM**, our framework for enhancing polygenic risk score (PGS) transferability through integration of genomic annotations, multi-ancestry modeling, and fine-mapping results.

The repo is organized into three main folders:

- **[`PRISM/`](./PRISM/)** — core code for running PRISM  
- **[`Analysis/`](./Analysis/)** — scripts for downstream analyses and figure generation  
- **[`Benchmarking/`](./Benchmarking/)** — implementations of alternative PGS approaches for comparison  

---

## 📂 Repository Structure

### 1. [`PRISM/`](./PRISM/)
Core methods for data curation, cV2F training, and running PRISM.

- **Installing PRISM**  
  - For background on basic model settings, see:  
    [*Power of inclusion: Enhancing polygenic prediction with admixed individuals*](https://www.sciencedirect.com/science/article/pii/S0002929723003518#sec3).
  - System Requirements: PRISM can be run on Linux, macOS, or Windows systems with a 64-bit processor
  - A complete software environment is provided in the yml file [PRISM_env.yml](https://github.com/lucy-tian/PRISM/blob/main/PRISM_env.yml)
  - External Binaries: The pipeline requires the following command-line tools to be available on your system PATH:
    - PLINK2
    - zstdcat
  ```bash
    # Download PRISM
    git clone https://github.com/lucy-tian/PRISM.git
    cd PRISM

    # Create and activate the PRISM environment
    conda env create -f prism_environment.yml
    conda activate prism_env

    # Add external binaries to PATH
    export PATH="/path/to/plink2_folder:$PATH"
    export PATH="/path/to/zstd_folder:$PATH"
  ```
  - Typical install time: 10–25 minutes on a standard desktop computer, depending mainly on R package compilation time.
    
- **Generating a penalty factor file**
  - PRISM requires an annotation table to generate the penalty factor assignment file.  
    The annotation table should look like:
      ```{r, eval=FALSE, indent="  "}
        #CHROM   POS     ID      REF     ALT     <annotations: cV2F_Blood  cV2F_Liver  ...>
        ...
      ```

    *Note:* The current script also includes additional annotation columns  
    **`Csq_group`**, **`ClinVar`**, **`geno_source`**, and **`in_LDSC_hm3`**  
    for extra penalty-factor assignment logic (see *Methods* section of the manuscript).  
    To modify penalty assignments for your own use case, edit the following section of the R script:

      ```r
      snpnet_w = case_when(
        (
          (Csq_group == "PTVs") |
          (ClinVar == "pathogenic")
        ) ~ 0.5,

        (
          (geno_source == "hla") |
          (Csq_group == "PAVs") |
          (ClinVar == "likely_pathogenic")
        ) ~ 0.6,

        (
          !!sym(paste('cV2F', arg1, sep = "_")) >= arg2
        ) ~ 0.7,

        (
          (!in_LDSC_hm3)
        ) ~ 1.2,

        TRUE ~ 1
      )
      ```

  - To generate the penalty factor file, run:
      ```bash
      Rscript snpnet_penalty_factors_tissue_specific.R <tissue_name (e.g., Blood)> <cV2F_score_cutoff> <path_to_annotation_file>
      ```
- **Running PRISM**  
  - PRISM requires specification of **9 flags**, which should be set in `PRISM/iPGS_pipeline/paths.sh`.

    **Three directory paths**

    - `data_d` – path to the local PRISM repository  
      ```bash
      data_d=/path/to/PRISM_folder
      ```
    - `snpnet_repo_d` – path to the  `snpnet` R package  
      ```bash
      snpnet_repo_d=/path/to/snpnet_repo
      ```
    - `p_factors` – directory where PRISM penalty files is stored  
      ```bash
      p_factors=/path/to/penalty_files_folder
      ```

    **Five input files**

    - `genotype_pfile` – prefix for genotype files in PLINK2 pfile format (`.pgen`, `.pvar`, `.psam`)  
      ```bash
      genotype_pfile=/path/to/genotype_prefix
      ```

    - `phenotype_file` – phenotype table with one row per individual  
      - Expected columns (tab-delimited):  
      ```{r, eval=FALSE, indent="  "}
        #FID  IID     population    split    <covariates: age, sex, ...>    <phenotype ID: INI1001, INI10002...>
        0    111     WB            train    <25, 0, ...>                   <163, 10.8, ...>
        0    112     NWB           test     <28, 1, ...>                   <178, 8.6, ...>
        0    113     Afr           split    <46, 1, ...>                   <183, 12.3, ...>
        ...
      ```

    - `pheno_names` – mapping from phenotype IDs to human-readable names  
      ```{r, eval=FALSE, indent="  "}
        #GBE_ID        GBE_NAME            N_available
        INI30120       Lymphocyte Count     406000
        INI50030700    eGFR                 398500
        ...

    - `keep_file` – `.keep` file listing individuals in the training + validation sets  

    - `sscore_afreq_ref` – allele-count reference file for the training population  
      - Typically generated from the training set (e.g. with `plink2 --freq counts`).  
      - Expected header (tab-delimited):
        ```{r, eval=FALSE, indent="  "}
        #CHROM  POS  ID  REF  ALT  HOM_REF_CT  HET_REF_ALT_CTS  TWO_ALT_GENO_CTS  HAP_REF_CT  HAP_ALT_CTS  MISSING_CT  OBS_CT
        ```

    **One covariate-specification string**

    - `covariates_type` – a single string listing covariate names separated by commas (no spaces).  
      ```bash
      covariates_type="age,sex,PC1,PC2,PC3,PC4,PC5"
      ```

    These 9 flags are read by the PRISM pipeline from `PRISM/iPGS_pipeline/paths.sh` and control where data are found and how models are fit.
    
  - To execute PRISM, run:
  ```bash
    cd iPGS_pipeline/
    sh PRISM_execute/fit_and_eval.sh <phenotype ID>
  ```
  
  - Expected run time: Expected run time: Typically 3–7 days on a normal desktop-class machine. Runtime depends heavily on R/SNPnet model fitting and available CPU/RAM.
    
- **PRISM Output**
  - PRISM output will be stored under the folder:  
    ```
    ${data_d}/PRISM_execute/fit_w_val/<phenotype_ID>/
    ```
    containing the following main files:

    - **`snpnet.BETAs.tsv.gz`**  
      SNP-level effect sizes for the fitted PRISM model.

    - **`snpnet.sscore.zst`**  
      Polygenic scores, stratified by `split` (train/val/test) and by `population`.

    - **`snpnet.eval.tsv.gz`**  
      PRISM model performance metrics (e.g., R²), evaluated across cohorts stratified  
      by population, dataset split, and other relevant groupings.

---

### 2. [`Analysis/`](./Analysis/)
Code for analyses and figures presented in the manuscript. Organized according to major sections:
- **Ancestry-specific analysis**  
- **Tissue-specific analysis**  
- **Biological interpretation**  
- **Benchmarking**  

Scripts here are mainly for **plot generation** and result visualization.

---

### 3. [`Benchmarking/`](./Benchmarking/)
Implementations of alternative PGS approaches used for comparison with PRISM.  

- **[IMPACT](https://github.com/immunogenomics/IMPACT/tree/master)**  
- **[SBayesRC](https://github.com/zhilizheng/SBayesRC)**  
- **[PolyPred](https://github.com/omerwe/polyfun)**  

Each subfolder contains scripts or wrappers for applying the corresponding method to our evaluation framework.

---

## Data & Resources

- **Data Curation**  
  - Details on variant-to-function (**cV2F**) scores and fine-mapping results:  
    [*A consensus variant-to-function score to functionally prioritize variants for disease*](https://www.biorxiv.org/content/10.1101/2024.11.07.622307v2.full.pdf).

- **Train cV2F**  
  - For replicating or retraining cV2F scores, see:  
    [Deylab999MSKCC/cv2f](https://github.com/Deylab999MSKCC/cv2f/tree/main).

- **Results**  
  - PGS performance for each PRISM model is reported in **Supplementary Table 8** of the manuscript.  
  - Derived SNP-level effect sizes for each PRISM model are hosted at Dropbox:  
    [SNP-level effect sizes](https://www.dropbox.com/scl/fo/9i3qzsolbaj79mn2bi4tn/AGHn_dQ6I6MdraGLk8dXZb4?rlkey=6j8cvbj21uy0wl0inkiwnp78o&st=nwvp40ud&dl=0).  
  - Data are organized by ancestry (*African*, *European*) and filenames follow the format:  
    ```
    [tissue]_[trait_id].tsv.gz
    ```
  - The `baseline/` folder contains results from a tissue- and ancestry-agnostic baseline model.  
  - Trait IDs correspond to phenotypes as follows:
    ```json
    {
      "INI30120"    : "Lymphocyte Count",
      "INI50030700" : "eGFR",
      "INI1003063"  : "FEV1 FVC_ratio",
      "INI20030780" : "LDL-C"
    }
    ```
## License

PRISM is released under the GNU Affero General Public License v3.0 (AGPL-3.0).

Under this license, any modifications or derivative works—including those deployed in networked or server environments—must also be released under the AGPL-3.0.

Organizations wishing to use PRISM privately or without these disclosure obligations may contact us to obtain a commercial license.

### Commercial Use

A commercial license is available for organizations that prefer to use PRISM
without AGPL-3.0 obligations. Please contact the corresponding authors of the relevant manuscript for details.


## 📖 References
- [*A consensus variant-to-function score to functionally prioritize variants for disease*](https://www.biorxiv.org/content/10.1101/2024.11.07.622307v2.full.pdf)  
- [*Power of inclusion: Enhancing polygenic prediction with admixed individuals*](https://www.sciencedirect.com/science/article/pii/S0002929723003518#sec3)  

---
