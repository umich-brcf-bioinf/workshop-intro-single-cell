This is the custom curriculum for the UM Bioinformatics Core Intro to Single Cell RNA-Seq Workshop.

* This workshop content is licensed under a [Creative Commons Attribution 4 License](https://creativecommons.org/licenses/by/4.0/).

* The workshop Code of Conduct has been adapted the NumFocus Code of Conduct (https://numfocus.org/code-of-conduct) which itself draws frin from numerous sources, including the Geek Feminism wiki, created by the Ada Initiative and other volunteers, which is under a Creative Commons Zero license, the Contributor Covenant version 1.2.0, the Bokeh Code of Conduct, the SciPy Code of Conduct, the Carpentries Code of Conduct, and the NeurIPS Code of Conduct.

* Links to web content:
  
  - [Main](https://umich-brcf-bioinf.github.io/workshop-intro-single-cell/main/html/)
  - [Most recent date branch](https://umich-brcf-bioinf.github.io/workshop-intro-single-cell/release/html/) 
  - [7/30/2025](https://umich-brcf-bioinf.github.io/workshop-intro-single-cell/2025-07-30/html/)
  - [4/16/2025](https://umich-brcf-bioinf.github.io/workshop-intro-single-cell/2025-04-16/html/)
  
* Input data for learners:
   - Copies of files pushed to learner's home directories are located in `/efs/workshop/isc/workshop_setup/shared_sample_data` on the AWS server

__For site knitting:__

* For repo setup, run `repo_setup.sh` script to ensure that untracked input files are copied to individual/local copy
	* Note that when working on LH/GL the untracked input files are pulled from the 'workshop' folder in ActiveProjects so those files need to check and updated to ensure they match the version(s) uploaded to AWS
* Before making changes, run `repo_clean.sh` to clear out the object caches and .RData files that are manually writen
	* This will remove `geo_so` files within `rdata` folder so the hidden code blocks are executed
	* Note that `cache` folders and `inputs` and `results` are not tracked as part of the github repo
* For editing and testing, `cd` into copy of repo and use singularity environment that's included as a comment at the top of the `repo_clean.sh`
* For first knit for branch, either:
	* Run `repo_setup.sh` script to pull (untracked) input files, then:
		* Execute first two lines of `knit_site.R` on command line (on Lighthouse/GL or AWS) OR
		* Open remote desktop and activate environment with rstudio 
* Run knitting of pages in order - a hidden starting block will read in Robject modified in the previous module
	* Note: table formatting has been updated to use kable styling so knitting will require installing both kable and klippy
