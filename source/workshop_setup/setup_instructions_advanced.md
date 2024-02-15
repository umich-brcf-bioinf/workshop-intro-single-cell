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


## Windows setup

### Installing R/RStudio (Windows)

1.  RStudio depends on the R programming environment, so we have to
    install that first. In a web browser, open:
    <https://cran.rstudio.com/bin/windows/base/>
    and click "Download R 4.3.2 for Windows" (the version may be
    slightly different). Open the downloaded executable to launch the R
    installer.

2.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes about 2 minutes; click **Finish** when prompted to close the
    installer.

3.  To install RStudio, in a web-browser, open:
    <https://rstudio.com/products/rstudio/download/#download>
    and click on **Download RStudio Desktop for Windows**. Open the
    downloaded executable to launch the installer.

4.  The installer will either prompt you to login as an Admin user or
    (if your current account has Admin privileges) simply ask you to
    allow it to make changes. Click **Yes**

5.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes less than one minute; click **Finish** when prompted to close
    the installer.

6.  Press Windows+R keys to open the **Run** dialog; type **"RStudio"**
    in the text box and hit enter. This will launch a new RStudio
    window. The RStudio window is divided into several panes. The lower
    left pane shows the **Console** tab and will show some text followed
    by a command prompt (\>):

> R version 4.3.2 (2023-10-31) -- "Eye Holes"
  Copyright (C) 2023 The R Foundation for Statistical Computing
  Platform: x86_64-pc-linux-gnu (64-bit)
  \
  R is free software and comes with ABSOLUTELY NO WARRANTY.
  You are welcome to redistribute it under certain conditions.
  Type 'license()' or 'licence()' for distribution details.
  \
  R is a collaborative project with many contributors.
  Type 'contributors()' for more information and
  'citation()' on how to cite R or R packages in publications.
  \
  Type 'demo()' for some demos, 'help()' for on-line help, or
  'help.start()' for an HTML browser interface to help.
  Type 'q()' to quit R.
  \


7.  The workshop exercises requires the installation of special R
    libraries. To install them into RStudio: 
    
    - one at a time, copy the blocks of code below and paste into the RStudio 
      **Console** tab
    - make sure the block you just pasted is highlighted and press **Enter** 
      to execute.
    - as the installation progresses, you might see red text flash by 
      and that's ok (typically an informative blurb or minor warning that has
      no downstream impact). 
      
```
if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {\
    install.packages(\"BiocManager\")\
    requireNamespace(\"BiocManager\", quietly = TRUE)
}
install.packages("devtools")
install.packages("remotes")
```

```
remotes::install_github("satijalab/seurat", "seurat5")
remotes::install_github("satijalab/seurat-data", "seurat5")
remotes::install_github("satijalab/seurat-wrappers", "seurat5")
```

```
devtools::install_github("immunogenomics/presto")
BiocManager::install("glmGamPoi")
```

8. Note: These installations automatically trigger the
    installation of a litany of dependent libraries so you will see
    repeated progress bars and code flying by in the Console window.
    The large set of dependencies mean that if you are installing from scratch
    **it can take several hours**. So now is a good time to make some
    coffee/tea (or a loaf of bread) while RStudio cooks.

9. If there was a problem during the installation, R will
    display a block like this:
    
> Warning in install.packages:
  Installation of package 'xxx' has non-zero exit status.

  If this comes up, you will need to troubleshoot the install based on the error
  or the specific package that failed.
  
10. If you don't see the above warning, and the output ends like below then you 
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

11. Press **Control-Q** close RStudio; when prompted to *Save workspace
    image...*, click **Don't Save**.

#### Your Windows workstation is ready for the workshop. Thank you for your patience and fortitude.

### Notes (Windows)

-   Following the workshop, you can remove R and
    RStudio. As an Admin user, go Start \> Settings \> Apps & Features.
    Click on the program to remove and click Uninstall.

## Macintosh setup

### Installing R/RStudio (Macintosh)

1.  RStudio depends on the R programming environment, so we have to
    install that first. In a web browser, open:

-   <https://cran.rstudio.com/bin/macosx/>

    and click the link "R-4.3.2.pkg" (the version may be slightly
    different). Open the downloaded executable to launch the R
    installer.

2.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Continue**) and when prompted,
    click **Install**. The installer will prompt you to confirm your
    username/password. The installer will show a progress bar; the
    process takes about 1 minutes; click **Finish** when prompted to
    close the installer.

3.  To install RStudio, in a web-browser, open:
-   <https://rstudio.com/products/rstudio/download/#download>
    and click on **Download RStudio Desktop for Mac**.

4.  Opening the downloaded executable opens a window with "RStudio" and
    your Applications folder icons. Drag the RStudio into the
    Applications folder. (If you see a dialog that claims RStudio
    already exists, click **Keep Both**.) The installer will prompt you
    to confirm your username/password.

5.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes less than one minute.

6.  When completed, open the Applications folder and double-click on the
    RStudio application. You will see a dialog
    "*RStudio.app is an app downloaded from Internet. Are you sure you
     want to open it?*"
    Click **Open**.
    This will launch a new RStudio window. The RStudio window is divided
    into several panes. The lower left pane shows the **Console** tab and
    will show some text followed by a command prompt (\>):

> R version 4.3.2 (2023-10-31) -- "Eye Holes"
  Copyright (C) 2023 The R Foundation for Statistical Computing
  Platform: x86_64-pc-linux-gnu (64-bit)
  \
  R is free software and comes with ABSOLUTELY NO WARRANTY.
  You are welcome to redistribute it under certain conditions.
  Type 'license()' or 'licence()' for distribution details.
  \
  R is a collaborative project with many contributors.
  Type 'contributors()' for more information and
  'citation()' on how to cite R or R packages in publications.
  \
  Type 'demo()' for some demos, 'help()' for on-line help, or
  'help.start()' for an HTML browser interface to help.
  Type 'q()' to quit R.
  \


7.  The workshop exercises requires the installation of special R
    libraries. To install them into RStudio: 
    
    - one at a time, copy the blocks of code below and paste into the RStudio 
      **Console** tab
    - make sure the block you just pasted is highlighted and press **Enter** 
      to execute.
    - as the installation progresses, you might see red text flash by 
      and that's ok (typically an informative blurb or minor warning that has
      no downstream impact). 
      
```
if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {\
    install.packages(\"BiocManager\")\
    requireNamespace(\"BiocManager\", quietly = TRUE)
}
install.packages("devtools")
install.packages("remotes")
```

```
remotes::install_github("satijalab/seurat", "seurat5")
remotes::install_github("satijalab/seurat-data", "seurat5")
remotes::install_github("satijalab/seurat-wrappers", "seurat5")
```

```
devtools::install_github("immunogenomics/presto")
BiocManager::install("glmGamPoi")
```

8. Note: These installations automatically trigger the
    installation of a litany of dependent libraries so you will see
    repeated progress bars and code flying by in the Console window.
    The large set of dependencies mean that if you are installing from scratch
    **it can take several hours**. So now is a good time to make some
    coffee/tea (or a loaf of bread) while RStudio cooks.

9. If there was a problem during the installation, R will
    display a block like this:
    
> Warning in install.packages:
  Installation of package 'xxx' has non-zero exit status.

  If this comes up, you will need to troubleshoot the install based on the error
  or the specific package that failed.
  
10. If you don't see the above warning, and the output ends like below then you 
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

11. Press **Control-Q** close RStudio; when prompted to *Save workspace
    image...*, click **Don't Save**.

####  Your Macintosh workstation is ready for the workshop. Thank you for your patience and fortitude.

### Notes (Macintosh)

-   Following the workshop, you can remove R and RStudio. Open
    your Applications directory, and drag the R and RStudio application
    into the Trash.
