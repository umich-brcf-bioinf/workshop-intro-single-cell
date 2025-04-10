---
title: "Introduction to Single-Cell Analysis"
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
<style type="text/css">
body{ /* Normal  */
      font-size: 14pt;
  }
</style>
# Check your Prereq

Thank you for your interest in the upcoming **Intro to Single-Cell Analysis** workshop. 


Learners in our previous workshop sessions have emphasized that **familiarity with R syntax and the RStudio interface** really helps you follow along and get the most out of the workshop.

---

## Do I need to brush up on R/R-Studio?

To confirm that you have the requisite background, you might consider the list of tasks below.

In an R environment:

- Create a data frame from a comma separated text file.
- Given a data frame with columns **city, state, zipcode**, create a smaller data frame with all rows but only columns **city** and **zipcode**.
- Given a data frame with columns **city, state, zipcode** create a smaller data frame with all columns but only rows where **state** is **Michigan**.
- Given a data frame with columns **city, state, zipcode** create a smaller data frame with a count of unique **zipcodes** for each **state**.


If you can provide the syntax for at least 3 of the above tasks off the top of your head, you’re in good shape. Otherwise, you may want to review key R/R-Studio techniques.

---

## How can I brush up?

- **Work through a self-guided tutorial:**
  - Checkout our [Intro to R and RStudio ](https://umich-brcf-bioinf.github.io/workshop-intro-r-rstudio/main/html/r-01-introduction.html){target="_blank"} workshop lesson.
  - Self-guided [R training from the Software Carpentry Institute](https://software-carpentry.org/lessons/){target="_blank"}
  - Peruse the first section of [R for Data Science](https://r4ds.hadley.nz/data-visualize){target="_blank"} (Hadley Wickham et al.).

  
- **Watch a video:**
  - [An introduction to the R programming language for Bioinformatics students](https://www.youtube.com/watch?v=bekFrlW0gww){target="_blank"}.<br/>
    (Sample data [here](https://drive.google.com/drive/folders/1mOCELXFb-b91C9mvfb2zD9nUTvNqqihO?usp=share_link){target="_blank"}.)

- Also, if you need a practice R/R-Studio environment, [PositCloud](https://posit.cloud/){target="_blank"} allows you to run R-Studio over the web for free. Click here to [login to PositCloud](https://posit.cloud/content/10156998){target="_blank"} and launch an R session.

---

Happy learning - and please let us know if you have questions!

[bioinformatics-workshops@umich.edu](mailto:bioinformatics-workshops@umich.edu)
