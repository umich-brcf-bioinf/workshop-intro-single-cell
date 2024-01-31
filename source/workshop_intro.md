---
title: "Intro to Single Cell RNA-Seq Workshop"
author: "UM Bioinformatics Core Workshop Team"
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
</style>

# Workshop goals and topics

 - By the end of the workshop, attendees will be able to
    - Describe (at a high level) how sc samples are sequenced along with strengths of a few popular library preps.
    - Create and interpret preliminary QC visualizations from a sc experiment.
    - Use Seurat, PCA, and UMAP to create cell-type clusters and projections.
    - Identify marker genes and annotate clusters based on gene expression.
    - Execute and visualize differential expression across clusters.


 - This workshop is targeted toward researchers who would like to be able to run
   scRNA-Seq analysis independently. This workshop assumes some experience working with 
   R/R-Studio.  Basic understanding of NGS sequencing protocols and/or familiarity with
   bulk RNA-Seq analysis is desirable.

 - Please let us know if there is anything we can do to improve the workshop experience.

 - Our purpose is not to be exhaustive, there is a lot that we cannot cover in
   the allotted time, and we don't expect anyone to be an expert at the end of
   the workshop. But we hope you will have a familiarity with key concepts, data types,
   tools, and how they all connect to one another in the service of a biological question.


## About the workshop team
| ![](images/headshots/headshot_cgates.jpg) | ![](images/headshots/headshot_mbradenb.jpg) | ![](images/headshots/headshot_rcavalca.jpg) |
|:-:|:-:|:-:|:-:|
| **Chris** | **Marci** | **Raymond** |
![](images/headshots/headshot_damki.jpg) | ![](images/headshots/headshot_ncarruth.jpg) | ![](images/headshots/headshot_trsaari.jpg) |
| **Dana** | **Nick** | **Travis** |
![](images/headshots/headshot_mkdohert.jpg) | ![](images/headshots/headshot_okoues.jpg) | ![](images/headshots/headshot_phillipa.jpg) |
| **Matt** | **Olivia** | **Tricia** |


<br/>

# Code of Conduct

- Be kind to others. Do not insult or put down others. Behave professionally. Remember that
  harassment and sexist, racist, or exclusionary jokes are not appropriate for the workshop.

- All communication should be appropriate for a professional audience including people of many
  different backgrounds. Sexual language and imagery is not appropriate.

- The Bioinformatics Core is dedicated to providing a harassment-free community for everyone,
  regardless of gender, sexual orientation, gender identity and expression, disability, physical
  appearance, body size, race, or religion. We do not tolerate harassment of participants in any
  form.

- Thank you for helping make this a welcoming, friendly community for all.

- If you have questions about the CoC please reach out to the hosts during the workshop, or
  email us at bioinformatics-workshops@umich.edu.

- To report a CoC incident/concern, please email Chris Gates (Bioinformatics Core, Managing
  Director) at cgates@umich.edu or contact the University of Michigan Office of Institutional
  Equity at institutional.equity@umich.edu.
<br/>
<br/>

# Using Zoom and Slack

- We will be recording this session. Recordings will be available to participants  
  following the workshop.


## <img src="images/workshop_intro/zoom_logo.png" alt="Zoom" width=120/>

- Zoom controls are at the bottom of the Zoom window:
![Zoom controls](images/workshop_intro/zoom_controls.png?s)

- To minimize distractions, we encourage participants to keep their audio muted
  (unless actively asking a question).
- To maximize engagement, we encourage participants to keep their video on.
- Slack works better than Zoom's Chat function so avoid Zoom Chat for now.
- You can enable transcription subtitles for your view.
- We will be using Breakout Rooms occasionally for ad-hoc 1-1 helper support.
  We will review this in detail together in a few minutes.
- Zoom's "non-verbal controls" are a useful way to interact
  - Depending on your version of Zoom you can access these either
    - in the **Reactions** button on you main Zoom window
    - at the bottom of the **Participant** pane

![Zoom non verbals](images/workshop_intro/zoom_nonverbals.png)

  - **Raise Hand** to request clarification or ask a question. (Same an in-person workshop.)
  - Instructors will use **Green check** and **Red X** to poll the group
    at checkpoints along the way.
<br/>
<br/>

### Exercise: Use Zoom non-verbals
  - Everyone use Zoom to raise your hand. <img src="images/workshop_intro/zoom_raise_hand.png" width=50/>
<br/>
<br/>


### Exercise: Using Zoom Breakout Rooms
Take a moment to briefly introduce yourself (name, dept/lab, area of study) in a
breakout room.

  - Zoom: Click Breakout Rooms
  - Find the room corresponding to the first letter of your last name
  - Hover over the number to the right and click **Join**.
  - When you have completed introductions, you can leave the breakout room to rejoin the main room.

  ![Breakout rooms](images/workshop_intro/zoom_breakout_room.png)
<br/>
<br/>


## <img src="images/workshop_intro/slack_logo.png" alt="Slack" width=120/>

- **Slack** can be used to communicate to the group or to individuals and has a
  few features/behaviors that we prefer over Zoom's Chat functionality.

![Posting messages to Slack](images/workshop_intro/slack.png)

  - Slack messages will be posted to the **2024-02-intro-single-cell** channel.  
    Click on the channel in the left pane (1) to select this channel.
  - You can type in the **message field** (2); click **send** (3) to
    post your message to everyone.
  - Helpers will respond in a Slack thread (or pose the question to the instructor)
  - You can respond in a message thread by hovering over a message to trigger
    the message menu and clicking the **speech bubble** (4).


### Exercise: Responding in Slack thread
In your scRNA-Seq experiment, what do/did you hope to find? 
<br/>
<br/>


## Review of Key communication patterns
| | <img src="images/workshop_intro/zoom_logo.png" alt="Zoom" width=120/> | <img src="images/workshop_intro/Slack_logo.png" alt="Slack" width=120/> |
|-|:-:|-|
| "I have an urgent question" | <img src="images/workshop_intro/zoom_raise_hand.png" width=50/> | **Post a question** |
| "I have a general question" | | **Post a question** |
| "I'm stuck / I need a hand" | | **Post a note** |
| Instructor check-in | <img src="images/workshop_intro/zoom_green_check.png" width=50/> -or- <img src="images/workshop_intro/zoom_red_x.png" width=50/> | |
| Instructor Slack question | | **Respond in Slack thread** |
<br/>
<br/>

### Exercise: Group checkpoint
  - Using Zoom, give me a **green-check** if you feel like you understand
    communication patterns or **red-X** if you need clarification.
<br/>
<br/>

## Arranging your screens
It is important that you can see:

 - Zoom (instructor's shared screen + reactions)
 - Your R/Studio screen
 - Slack
 - Lesson plan web page

![arranged screens](images/workshop_intro/arranged_screens.png?)
<br/>
<br/>



### Workshop topics at a glance
![](images/wayfinder/wayfinder_workshop_intro.png)

# Any questions?
 -
 -


---

# Thank you to our sponsors/contributors

![](images/workshop_intro/sponsor_logos.png)

## [UM BRCF Bioinformatics Core](https://brcf.medicine.umich.edu/bioinformatics){target="_blank"}

- The University of Michigan BRCF Bioinformatics Core is a team of analysts that help
  researchers **design, analyze, and interpret** high-throughput genomics experiments.
- Last year we helped about 60 researchers design and execute about 100 projects
  including gene expression, epigenetic, variant identification, functional
  enrichment and many other kinds of analyses.
- We provide letters of support for grant proposals.
- We are creating a series of bioinformatics-focused workshops.

## [University of Michigan Library](https://www.lib.umich.edu/){target="_blank"}

  Our mission is to support, enhance, and collaborate in the instructional, research, and
  service activities of faculty, students, and staff, and contribute to the common good by
  collecting, organizing, preserving, communicating, sharing, and creating the record of human
  knowledge.

## [Biomedical Research Core Facilities](https://brcf.medicine.umich.edu/){target="_blank"}

  Biomedical Research Core Facilities (BRCF) helps researchers economically take advantage of
  the latest technology and collaborate with top experts in the field. Established in 1986, the
  BRCF was formed to offer centralized access to research services and equipment.
<br/>
<br/>


## Acknowledgements

* TODO Acknowledgements

* The workshop Code of Conduct has been adapted the NumFocus Code of Conduct (https://numfocus.org/code-of-conduct) which itself draws frin from numerous sources, including the Geek Feminism wiki, created by the Ada Initiative and other volunteers, which is under a Creative Commons Zero license, the Contributor Covenant version 1.2.0, the Bokeh Code of Conduct, the SciPy Code of Conduct, the Carpentries Code of Conduct, and the NeurIPS Code of Conduct.

* This workshop content is licensed under a [Creative Commons Attribution 4 License](https://creativecommons.org/licenses/by/4.0/).
