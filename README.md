This is the custom curriculum for the UM Bioinformatics Core Intro to Single Cell RNA-Seq Workshop. 


* This workshop content is licensed under a [Creative Commons Attribution 4 License](https://creativecommons.org/licenses/by/4.0/).

* The workshop Code of Conduct has been adapted the NumFocus Code of Conduct (https://numfocus.org/code-of-conduct) which itself draws frin from numerous sources, including the Geek Feminism wiki, created by the Ada Initiative and other volunteers, which is under a Creative Commons Zero license, the Contributor Covenant version 1.2.0, the Bokeh Code of Conduct, the SciPy Code of Conduct, the Carpentries Code of Conduct, and the NeurIPS Code of Conduct.

__For site knitting:___

* For first knit for branch, either:
	* Run `repo_setup.sh` script to pull (untracked) input files, then:
		* Execute first two lines of `knit_site.R` on command line (on comps/Lighthouse or AWS) OR
		* Open remote desktop and activate environment with rstudio 
* Run knitting of pages in order - a hidden starting block will read in Robject modified in the previous module
* After modifying an individual page, `geo_so` must be removed before re-knitting so the hidden code block is executed
	* The `repo_clean.sh` script can be used to clear all inputs and reset state 