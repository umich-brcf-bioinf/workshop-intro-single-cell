---
title: "Intro to scRNA-Seq Workshop: Advanced setup instructions"
author: "UM Bioinformatics Core Workshop Team"
output:
        html_document:
            theme: paper
            toc: true
            toc_depth: 6
            toc_float: true
            number_sections: false
            fig_caption: false
            markdown: GFM
            code_download: false
---

- The steps below are oriented toward advanced users who would like to install
  and configure select software on their local workstation. These steps are not
  necessary or recommended, but may provide context for advanced use-cases.

- These instructions build on the install detailed in the [basic setup
  instructions](setup_instructions.html); please see that document for more
  details on setting up Zoom, Slack, and also for getting help.

- Note that if you do not have Administrative privileges it will be tricky
  to install R/RStudio; you may need to coordinate with your
  System Admin/IT Support team to get these installed.

- Note that these instructions are focused on R/RStudio and dependent libraries 
  (e.g. Seurat). Seperate instructions on downloading and installing Cell Ranger
  can be found here: 
  https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in 

  #### Table of Contents
  -   [Windows setup](#windows-setup)
      -   Installing R/RStudio
      -   Notes
  -   [Macintosh setup](#macintosh-setup)
      -   Installing R/RStudio
      -   Notes
  - Installing required libraries
      - Conda & Conda packages
      - R packages
      

## Windows setup

### Installing R/RStudio (Windows)

1.1.  RStudio depends on the R programming environment, so we have to
    install that first. In a web browser, open:
    <https://cran.rstudio.com/bin/windows/base/>
    and click "Download R 4.3.2 for Windows" (the version may be
    slightly different). Open the downloaded executable to launch the R
    installer.

1.2.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes about 2 minutes; click **Finish** when prompted to close the
    installer.

1.3.  To install RStudio, in a web-browser, open:
    <https://rstudio.com/products/rstudio/download/#download>
    and click on **Download RStudio Desktop for Windows**. Open the
    downloaded executable to launch the installer.

1.4.  The installer will either prompt you to login as an Admin user or
    (if your current account has Admin privileges) simply ask you to
    allow it to make changes. Click **Yes**

1.5.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes less than one minute; click **Finish** when prompted to close
    the installer.

1.6.  Press Windows+R keys to open the **Run** dialog; type **"RStudio"**
    in the text box and hit enter. This will launch a new RStudio
    window. The RStudio window is divided into several panes. The lower
    left pane shows the **Console** tab and will show some text followed
    by a command prompt (\>):

```
R version 4.4.2 (2024-10-31 ucrt) -- "Pile of Leaves"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 

```

1.7.  Following the workshop, you can remove R and
    RStudio. As an Admin user, go Start \> Settings \> Apps & Features.
    Click on the program to remove and click Uninstall.


## Macintosh setup

### Installing R/RStudio (Macintosh)

1.1.  RStudio depends on the R programming environment, so we have to
    install that first. In a web browser, open: 
    [https://cran.rstudio.com/bin/macosx](https://cran.rstudio.com/bin/macosx){target="_blank"}. 
    - Note that newer Macs will use Silicon build, older Macs (pre 2020) will use Intel builds. 
      (You can confirm your architecture under About this Mac: Chip; "Apple M" indicates Silicon. )
    - Right-click to download the latest version (e.g. "R-4.4.2"", though the 
    specific version may be slightly different). 
    - Open the downloaded executable to launch the R installer.

1.2.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Continue**) and when prompted,
    click **Install**. The installer will prompt you to confirm your
    username/password. The installer will show a progress bar; the
    process takes about 1 minutes; click **Finish** when prompted to
    close the installer.

1.3.  To install RStudio, in a web-browser, open:
-   <https://rstudio.com/products/rstudio/download/#download>
    and click on **Download RStudio Desktop for Mac**.

1.4.  Opening the downloaded executable opens a window with "RStudio" and
    your Applications folder icons. Drag the RStudio into the
    Applications folder. (If you see a dialog that claims RStudio
    already exists, click **Keep Both**.) The installer will prompt you
    to confirm your username/password.

1.5.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes less than one minute.

1.6.  When completed, open the Applications folder and double-click on the
    RStudio application. You will see a dialog
    "*RStudio.app is an app downloaded from Internet. Are you sure you
     want to open it?*"
    Click **Open**.
    This will launch a new RStudio window. The RStudio window is divided
    into several panes. The lower left pane shows the **Console** tab and
    will show some text followed by a command prompt (\>):


```
R version 4.4.2 (2024-10-31) -- "Pile of Leaves"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin20

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

    Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 

```

1.7.  Following the workshop, you can remove R and RStudio. Open
    your Applications directory, and drag the R and RStudio application
    into the Trash.


### Installing required libraries

#### Install software dependencies

A few R libraries depend on external software programs.

2.1 Install the Miniconda package manager as detailed [here](https://docs.anaconda.com/miniconda/install/){target="blank"}

2.2 Open a command/teminal window and create and activate a new Conda environment 
```
conda create -n isc
conda activate isc
```

2.3 Install required software using Conda:

```
# Required by BPCells
conda install -y conda-forge::hdf5 

# Optional: used by Leiden clustering algorithm
conda install -y conda-forge::leidenalg

# Optional: used to convert Illumina BCLs to FASTQs
conda install -y bioconda::bcl2fastq-nextseq

```

2.4 Optional: Install Cell Ranger (Useful only if you have to rerun Cell Ranger on data you received from a sequencing facility.)
See detailed [installation instructions](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in){target="blank"} at 10x Genomics 


#### Install R packages
3.1  The workshop exercises requires the installation of special R
    packages. To install them into RStudio: 
    
    - one at a time, copy the blocks of code below and paste into the RStudio 
      **Console** tab
    - make sure the block you just pasted is highlighted and press **Enter** 
      to execute.
    - as the installation progresses, you might see red text flash by 
      and that's ok (typically an informative blurb or minor warning that has
      no downstream impact). 
      
```
# -------------------------------------------------------------------------
install.packages("devtools")
install.packages("remotes")

if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {\
    install.packages(\"BiocManager\")\
    requireNamespace(\"BiocManager\", quietly = TRUE)
}

# -------------------------------------------------------------------------
remotes::install_github("satijalab/seurat", "seurat5")
remotes::install_github("satijalab/seurat-data", "seurat5")
remotes::install_github("satijalab/seurat-wrappers", "seurat5")

# -------------------------------------------------------------------------
# Allows data to live on disk instead of in RAM
# https://github.com/bnprks/BPCells
remotes::install_github("bnprks/BPCells/r")


# -------------------------------------------------------------------------
# Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data
# https://github.com/ZJUFanLab/scCATCH
install.packages("scCATCH")

# -------------------------------------------------------------------------
# Differential gene expresion testing
# https://bioconductor.org/packages/DESeq2/
BiocManager::install("DESeq2", update = FALSE)

# =========================================================================
# OPTIONAL libraries
# Some of the libraries make the code execute faster and some enable some 
# off-menu options (e.g. Leiden clustering). The workshop notes will work fine 
# without these libraries below but may be slower.
# =========================================================================

# -------------------------------------------------------------------------
# Presto makes it fast and easy to run Wilcoxon rank sum test and auROC analysis on large datasets
# https://github.com/immunogenomics/presto
devtools::install_github("immunogenomics/presto")

# -------------------------------------------------------------------------
# Accelerates DESeq2 diffex analysis by fiting Gamma-Poisson Generalized Linear Models Reliably
# https://github.com/const-ae/glmGamPoi
BiocManager::install("glmGamPoi")

# -------------------------------------------------------------------------
# An alternative to Louvain clustering algorithm
# https://cran.r-project.org/package=leiden ; https://github.com/TomKellyGenetics/leiden
install.packages("leiden")
```

3.2 Note: These installations automatically trigger the
    installation of a litany of dependent libraries so you will see
    repeated progress bars and code flying by in the Console window.
    The large set of dependencies mean that if you are installing from scratch
    **it can take several hours**. So now is a good time to make some
    coffee/tea (or a loaf of bread) while RStudio cooks.

3.3. If there was a problem during the installation, R will
    display a block like this:
    
> Warning in install.packages:
  Installation of package 'xxx' has non-zero exit status.

  If this comes up, you will need to troubleshoot the install based on the error
  or the specific package that failed.
  
3.4. If you don't see the above warning, and the output ends like below then you 
   the libraries were successfully installed:

> The downloaded binary packages are in
  C:\\Users\\some\\path\\Temp\\downloaded_packages

   You can check this with the following code block:
```
# Quietly load Seurat and print the version
# It should be version 5 or later
suppressPackageStartupMessages(library(Seurat))
sessionInfo()$otherPkgs$Seurat$Version
```   

3.5. Press **Control-Q** close RStudio; when prompted to *Save workspace
    image...*, click **Don't Save**.

#### Your workstation is ready for the workshop. Thank you for your patience and fortitude.
