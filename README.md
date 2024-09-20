This is the custom curriculum for the UM Bioinformatics Core Intro to Single Cell RNA-Seq Workshop. 


* This workshop content is licensed under a [Creative Commons Attribution 4 License](https://creativecommons.org/licenses/by/4.0/).

* The workshop Code of Conduct has been adapted the NumFocus Code of Conduct (https://numfocus.org/code-of-conduct) which itself draws frin from numerous sources, including the Geek Feminism wiki, created by the Ada Initiative and other volunteers, which is under a Creative Commons Zero license, the Contributor Covenant version 1.2.0, the Bokeh Code of Conduct, the SciPy Code of Conduct, the Carpentries Code of Conduct, and the NeurIPS Code of Conduct.

__For site knitting:___

* For repo setup, run `repo_setup.sh` script to ensure that untracked input files are copied to individual/local copy
* Before making changes, run `repo_clean.sh` to clear out the object caches and .RData files that are manually writen
	* This will remove `geo_so` files within `rdata` folder so the hidden code blocks are executed
	* Note that `cache` folders and `inputs` and `results` are not tracked as part of the github repo
* For editing and testing, `cd` into copy of repo and use singularity environment that's included as a comment at the top of the `repo_clean.sh`
* For first knit for branch, either:
	* Run `repo_setup.sh` script to pull (untracked) input files, then:
		* Execute first two lines of `knit_site.R` on command line (on comps/Lighthouse or AWS) OR
		* Open remote desktop and activate environment with rstudio 
* Run knitting of pages in order - a hidden starting block will read in Robject modified in the previous module
