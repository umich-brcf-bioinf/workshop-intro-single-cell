## Setup input directory, assuming you're at repo root
mkdir -p source/inputs

## Pull lastest version of inputs (with updated CellRanger) - 1.28.25
## Note that path to input data will be different if we're knitting on the AWS server
cp -r /nfs/mm-isilon/bioinfcore/ActiveProjects/workshop/inputs/* source/inputs 
# -l flag triggered failure for updated input data 
