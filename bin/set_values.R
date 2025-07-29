workshop_vars = list(
  ## Common
  rstudio_server_url = "http://bfx-workshop02.med.umich.edu",
  ## Code of Conduct
  coc_contact = "University of Michigan Equity, Civil Rights, and Title IX Office",
  coc_contact_email = "ecrtoffice@umich.edu",
  ## Slack
  slack_channel = "2025-07-intro-single-cell",
  ## Wrap-up
  rstudio_server_enddate = "8/15/2025",
  ssh_download_dns = "bfx-workshop01.med.umich.edu",
  aws_s3_bucket = "https://umich-brcf-bioinf-workshop.s3.us-east-1.amazonaws.com",
  aws_s3_file = "ISC/workshop_isc_inputs-20250730.tgz",
  satija_scgd = "https://satijalab.org/scgd25/"
)

# The sessionInfo outputs will only be refreshed if this var exists and we 
# are running on AWS
on_aws = dir.exists('/efs/workshop')

