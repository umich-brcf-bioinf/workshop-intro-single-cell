---
title: "Introduction to Single-Cell Analysis: Setup instructions"
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

This document guides you through the installation of software necessary
to participate in the workshop. To complete
this setup, you will need:

-   A Macintosh or Windows workstation connected to the internet.
-   An email invitation to the workshop Slack workspace.
-   An email titled **"UM BioinfCore Workshop Login"** containing your
    individual login and password for the shared workshop server
-   About 10 minutes.


## Introduction

-   The virtual workshop will be conducted using Zoom and Slack for all days.
-   Attendees will use their own workstation/laptop to participate; the
    workstation/laptop should have a microphone, camera, and reliable
    internet access. While not required, for laptop users we recommend using a
    larger monitor for the optimal workshop experience.
-   The workshop will use Slack for group messaging. Messages posted to the
    general channel will be visible to all participants in the workshop.
-   For simplicity, participants will log in to a pre-built server which
    has the necessary data, software, and configuration already installed.
-   While not required or recommended, advanced users may elect to install some software
    locally; see [advanced setup instructions](setup_instructions_advanced.html).

## How to get help

While we have endeavored to make this setup process work robust and
comprehensive, installing bioinformatics software is tricky and we would
be happy to lend a hand to get things working.

-   If you have problems/questions, please don't hesitate to email us
    at:
    [bioinformatics-workshops@umich.edu](mailto:bioinformatics-workshops@umich.edu)

-   When emailing it will speed things along if you could include:
    -   Whether you are using Windows or Mac (and optionally which
        version of the OS you are using).
    -   Whether you have Administrative privileges on your workstation.
    -   The specific text of any error messages, if applicable.

## Setup steps

### **Zoom**

   If you have not used Zoom before, please use the following link to
   install "Zoom client for Meetings"\
   https://zoom.us/download

   - Additional Zoom settings for Macintosh:

     1.1  To enable screen sharing (useful for breakout rooms and tech support)

        - System Preferences \>\> Security & Privacy: click on the Privacy tab
        - Select **Screen Recording** on left tab
        - Scroll to the bottom of the right tab and verify **Zoom** is checked.

     1.2.  To enable remote control (useful for breakout rooms and tech support)

        - System Preferences \>\> Security & Privacy: click on the Privacy tab
        - Select **Accessibility** on left tab
        - Scroll to the bottom on the right tab and verify **Zoom** is checked.
        - If it's not checked, click the lower left lock icon and enter
            user password when prompted. You can now check Zoom in the
            right tab.
        - Click the lock again.

      1.3.  Close and relaunch Zoom client.

### **Slack**

   2.1. In the Slack email invitation from the hosts, click the link to join the Slack workspace.

   2.2. Follow the Slack instructions to login or create an account. You may
      edit your profile to add a photo if you choose.

   2.3. Note that once your Slack account is established and you have accepted the
      Slack invite you can access the Slack channel by pointing your web browser
      here: https://umbioinfcoreworkshops.slack.com

### **Launch R/RStudio (from the web)**

  3.1 Open the following web address: [https://bfx-workshop02.med.umich.edu](https://bfx-workshop02.med.umich.edu){target="_blank"}

  3.2 When prompted, enter the  username and password supplied by your hosts and
   click **Sign In**.

  3.3 The login page will immediately go blank - **it may take several seconds for the next page to appear** 
   (sometimes up to **60 seconds**, but usually faster). After this delay, you 
   should see a new webpage that looks something like this:

   ![RStudio initial screen](images/setup_instructions_rstudio.png)

   If you see a screen like that above, you have successfully logged into RStudio.

## Congratulations
Your workstation and logins are ready for the workshop; thank you for
   taking time to set this up in advance.

   - You can close your RStudio session by clicking the small orange "power"
     button in the top-right corner (highlighted in the picture above by a red
     arrow).  The webpage will announce **R Session Ended**; you may now close
     the web page.

### A note on the shared compute environment

Please note that this environment is optimized for the exercises in this
particular workshop but is likely unsuitable for analyzing your own
datasets. In particular:
- It is not sized for compute intensive operations or large storage.
- It is not secured for sensitive data of any kind.
- This environment is temporary and will be removed shortly after the
  conclusion of the workshop.

If you had problems or questions about the installation steps
please see the **[How to get help](#how-to-get-help)** section above for
more assistance.
