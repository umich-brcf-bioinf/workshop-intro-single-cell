# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("biomaRt","DESeq2"), update=FALSE, ask=FALSE)
# #
# missing <- setdiff(c("tidyr", "ggplot2", "pheatmap", "ggrepel", "formattable", "RColorBrewer", "matrixStats", "dplyr", "biomaRt", "DESeq2"), rownames(installed.packages()))
#
# if (!length(missing)) {
#   cat("Ready for Computational Foundations workshop\n")
#   } else {
#     cat("PROBLEM: could not install:", missing, "\n")
#   }
# install.packages("pheatmap")
# install.packages("ggrepel")
# install.packages("formattable")
# install.packages("tidyr")
# install.packages("RColorBrewer")
# install.packages("matrixStats")
# install.packages("dplyr")
# }

library(rmarkdown)

# The html from the files below don't have the nav bar

# TODO: review/build this content
render('source/workshop_setup/preregistration_info.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/preworkshop_checklist.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions.md', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions_advanced.md', output_dir='html/workshop_setup/')

# The html from the files below do have the nav bar, so if you make changes 
# that impact the navbar (e.g. file name changes) you should re-knit all of them.

render_site('source/index.md')
render_site('source/workshop_intro.md')

## add content pages
render_site('source/00A-OrientingOnScRNASeq.Rmd')
render_site('source/00B-CellRangerInAction.md')

render_site('source/01-GettingStarted.Rmd')
render_site('source/02-QCandFiltering.Rmd')
render_site('source/03-Normalization.Rmd')
render_site('source/04-PCAandIntegration.Rmd')
render_site('source/05-ProjectionAndClustering.Rmd')
render_site('source/06-MarkerVisualization.Rmd')
render_site('source/07-CellTypeAnnos.Rmd')
render_site('source/08-DifferentialExpression.Rmd')
render_site('source/09-IndependentExercise.Rmd')
render_site('source/00-ResourcesAndExtendedContent.Rmd')

render_site('source/workshop_wrap_up.md')
#clean_site(preview=TRUE)
