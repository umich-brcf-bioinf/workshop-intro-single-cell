## if not running on AWS - setup environment
# module load singularity
# singularity exec /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.11.0.sif R ## activate R
# singularity exec /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.11.0.sif rstudio ##activate Rstudio

## to ensure correct formatting for knit html, install required packages (klippy and kable)
# install.packages('remotes'); library(remotes)
# #remove.packages('klippy');
# remotes::install_github("umich-brcf-bioinf/workshop-klippy"); 
# 
#install.packages("kable")
#install.packages('devtools'); library(devtools)
#devtools::install_github("haozhu233/kableExtra")

library(rmarkdown)
library(klippy)
library(kableExtra)

# The html from the files below don't have the nav bar
render('source/workshop_setup/preregistration_info.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/preworkshop_checklist.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions.Rmd', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions_advanced.Rmd', output_dir='html/workshop_setup/')
render('source/workshop_setup/prereq_check.md', output_dir='html/workshop_setup/')


# The html from the files below do have the nav bar, so if you make changes 
# that impact the navbar (e.g. file name changes or reordering) you should 
# re-knit all of them.

render_site('source/index.md');
render_site('source/workshop_intro.Rmd');
render_site('source/instructor_cheatsheet.Rmd');
render_site('source/analysis_scripts.Rmd');

## content pages
render_site('source/00A-OrientingOnScRNASeq.Rmd');
render_site('source/01-GettingStarted.Rmd');
render_site('source/00B-CellRangerInAction.Rmd');
render_site('source/02-QCandFiltering.Rmd');
render_site('source/03-Normalization.Rmd');

render_site('source/04-PCAandIntegration.Rmd');
render_site('source/05-ProjectionAndClustering.Rmd');
render_site('source/clusters_faq.Rmd');

render_site('source/06-MarkerVisualization.Rmd');
render_site('source/07-CellTypeAnnos.Rmd');
render_site('source/08-DifferentialExpression.Rmd');
render_site('source/08A-AnalysisFinale.Rmd');

render_site('source/00-ResourcesAndExtendedContent.Rmd');
render_site('source/exercises.Rmd');

render_site('source/seurat-on-great-lakes.Rmd');
render_site('source/workshop_wrap_up.Rmd');
rm(list=ls(all.names = TRUE))
gc(verbose=TRUE, full=TRUE)
#clean_site(preview=TRUE)

