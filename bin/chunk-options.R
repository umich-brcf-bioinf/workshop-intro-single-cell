# These settings control the behavior of all chunks in the novice R materials.
# For example, to generate the lessons with all the output hidden, simply change
# `results` from "markup" to "hide".
# For more information on available chunk options, see
# http://yihui.name/knitr/options#chunk_options

library("knitr")

fix_fig_path <- function(pth) file.path(pth)


## We set the path for the figures globally below, so if we want to
## customize it for individual episodes, we can append a prefix to the
## global path. For instance, if we call knitr_fig_path("01-") in the
## first episode of the lesson, it will generate the figures in
## `fig/rmd-01-`
knitr_fig_path <- function(prefix) {
    new_path <- paste0('images/curriculum/',
                      prefix)
    opts_chunk$set(fig.path = new_path)
}

## We use the rmd- prefix for the figures generated by the lessons so
## they can be easily identified and deleted by `make clean-rmd`.  The
## working directory when the lessons are generated is the root so the
## figures need to be saved in fig/, but when the site is generated,
## the episodes will be one level down. We fix the path using the
## `fig.process` option.

opts_chunk$set(tidy = FALSE, results = "markup", comment = NA,
               fig.align = "center", fig.path = "images/curriculum/rmd-",
               fig.process = fix_fig_path,
               fig.width = 8.5, fig.height = 8.5,
               fig.retina = 2)
