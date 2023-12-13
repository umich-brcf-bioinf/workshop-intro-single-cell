---
title: "Intro to Single Cell RNA-Seq Workshop"
author: "UM Bioinformatics Core"
output:
        html_document:
            includes:
                in_header: header.html
            theme: paper
            toc: true
            toc_depth: 4
            toc_float: true
            number_sections: false
            fig_caption: true
            markdown: GFM
            code_download: false
---

<style type="text/css">

body, td {
   font-size: 18px;
}
code.r{
  font-size: 12px;
}
pre {
  font-size: 12px
}

a.external {
    background: url(images/external-link.png) center right no-repeat;
    padding-right: 13px;
}
</style>

## Wrapping up

We hope you now have more familiarity with key concepts, data types, tools, and how they all
connect to enable single-cell gene expression analysis from single cell RNA-Seq data.

![](images/wayfinder/wayfinder_workshop_wrap_up.png)

---

## Housekeeping

- TODO: Please take our optional <a class="external" href="https://forms.gle/rTDtL5ExtvhGryVs6" target="_blank">post-workshop survey</a> (5-10 minutes)

- TODO: 11/27 11-12 ET: <a class="external" href="https://forms.gle/Rfx1TQ64gw2TnizB9" target="_blank">Virtual office hours</a>

- We will email you a link to the final session recordings by next week.

- The website/notes for this workshop will be available.

- The <a class="external" href="https://umbioinfcoreworkshops.slack.com" target="_blank">UM Bioinformatics Core Workshop Slack channel</a> will be available
for 90 days.

---

## Looking ahead

#### Workshop environment
- TODO: RStudio workshop compute environment will be available until 11/30.
  - _Please save all your R scripts now_ so that we can "right-size" the compute environment immediately following today's workshop session.
- You can download files from the workshop environment from your terminal/command line window as below. 
  (You will need to substitute your actual workshop username and type workshop password when prompted.)
  ```
  mkdir intro_scrnaseq-workshop
  cd intro_scrnaseq-workshop
  scp -r YOUR_USERNAME@bfx-workshop01.med.umich.edu:"ISC*" .
  ```
  - TODO: Note that the full download is about 2-4Gb, so depending on your internet 
    speeds it could take while.

#### TODO Installing software locally
  - You can install necessary programs to run programs locally. Note
    that for typical data, the compute intensive steps (Day 1) assume your 
    computer has powerful compute (many CPUs and lots of RAM) and sizable 
    storage capacity. (i.e. It may not be practical to run these on your laptop.)
  - Installing bioinformatics software is non-trivial and comprehensive instructions
    to setup a complete compute environment are outside the scope of this workshop.
    (We are considering offering a Computational Reproducibility workshop that
    would cover this.) For the intrepid, see relevant links below:
    - [Installing Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)
    - [Conda environment YAML](https://github.com/umich-brcf-bioinf/2022-08-14-umich-rnaseq-demystified/blob/main/server_setup/conda_setup/workshop-rsd.yaml)
    - [Setting up R/R-Studio](workshop_setup/setup_instructions_advanced.html)


## University of Michigan Resources
- <a class="external"  href="https://datascience.isr.umich.edu/events/coderspaces/" target="_blank">UM CoderSpaces "office hours"</a> and UM CoderSpaces Slack workspace. _(See "Useful Resources" section in above page for instructions on how to join and access the CoderSpaces Slack workspace.)_
- Upcoming <a class="external" href="https://arc.umich.edu/events/" target="_blank">UM Advanced Research Computing workshops</a>.
- Advanced Research Computing (ARC) at University of Michigan hosts a
  high-performance computing (HPC) platform called _Great Lakes_ which combines
  high-end computers, fast/resilient storage, and pre-installed software.
  GreatLakes may be a good resource for folks who need to run the more compute
  intensive steps and a substantial block of compute and storage is subsidized
  by ARC making it essentially free to many researchers.
  - About <a class="external" href="https://arc.umich.edu/greatlakes/" target="_blank">Great Lakes HPC</a>.
  - About the ARC <a class="external" href="https://arc.umich.edu/umrcp/" target="_blank">Research Computing Package</a>.
  - Videos on <a class="external" href="https://www.mivideo.it.umich.edu/channel/ARC-TS%2BTraining/181860561/" target="_blank">getting started with Great Lakes</a>. (available to UM folks)


#### Training and support
  - Learning bioinformatic analyses more like a process than a task. Resources 
    to consider:
    - UM Bioinformatics Core links to <a href="https://brcf.medicine.umich.edu/cores/bioinformatics-core/training/" target="_blank">training resources/workshops</a>
    - TODO: Notes on Bash and R from the recent <a href="https://umich-brcf-bioinf.github.io/2023-07-31-umich-computational-foundations/html/" target="_blank">Computational Foundations Workshop</a>.<br/>
    - For more intro lessons and workshops in Bash / Git / R / Python : <a href="https://software-carpentry.org/lessons/" target="_blank">Software Carpentry</a> and the <a href="https://umcarpentries.org/" target="_blank">UM Software Carpentry Group</a>.
    - For more info on NGS analysis, we can highly recommend training materials from the <a href="https://hbctraining.github.io/main/" target="_blank">Harvard Chan Bioinformatics Core</a>.

---

## Thank you

![Sponsors](images/workshop_intro/sponsor_logos.png)

## Thank you to/from the workshop team
| ![](images/headshots/headshot_cgates.jpg) | ![](images/headshots/headshot_mbradenb.jpg) | ![](images/headshots/headshot_rcavalca.jpg) |
|:-:|:-:|:-:|:-:|
| **Chris** | **Marci** | **Raymond** |
![](images/headshots/headshot_damki.jpg) | ![](images/headshots/headshot_ncarruth.jpg) | ![](images/headshots/headshot_trsaari.jpg) |
| **Dana** | **Nick** | **Travis** |

<br/>



Thank you for participating in our workshop. We welcome your questions and feedback now and in the future.

Bioinformatics Workshop Team

[bioinformatics-workshops@umich.edu](mailto:bioinformatics-workshops@umich.edu) <br/>
[https://brcf.medicine.umich.edu/bioinformatics](https://brcf.medicine.umich.edu/bioinformatics)
