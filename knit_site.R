# module load singularity
# singularity exec /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.11.0.sif R

# install.packages('remotes'); library(remotes)
# remove.packages('klippy');remotes::install_github("umich-brcf-bioinf/workshop-klippy"); library(klippy)
library(rmarkdown)

# The html from the files below don't have the nav bar

render('source/workshop_setup/preregistration_info.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/preworkshop_checklist.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions_advanced.md', output_dir='html/workshop_setup/')

# The html from the files below do have the nav bar, so if you make changes 
# that impact the navbar (e.g. file name changes or reordering) you should 
# re-knit all of them.

render_site('source/index.md')
render_site('source/workshop_intro.md')


## add content pages
render_site('source/00A-OrientingOnScRNASeq.Rmd')
render_site('source/01-GettingStarted.Rmd')
render_site('source/00B-CellRangerInAction.md')
render_site('source/02-QCandFiltering.Rmd')
render_site('source/03-Normalization.Rmd')
render_site('source/04-PCAandIntegration.Rmd')
render_site('source/05-ProjectionAndClustering.Rmd')
render_site('source/06-MarkerVisualization.Rmd')
render_site('source/07-CellTypeAnnos.Rmd')
render_site('source/08-DifferentialExpression.Rmd')
render_site('source/09-IndependentExercise.Rmd')
render_site('source/00-ResourcesAndExtendedContent.Rmd')

render_site('source/workshop_wrap_up.Rmd')
rm(list=ls())
gc(verbose=TRUE, full=TRUE)
#clean_site(preview=TRUE)





