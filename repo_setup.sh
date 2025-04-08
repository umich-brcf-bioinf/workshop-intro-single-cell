## Setup input directory, assuming you're at repo root
mkdir -p source/inputs/prepared_data

## Pull lastest version of inputs (with updated CellRanger) - 1.28.25
## Note that path to input data will be different if we're knitting on the AWS server
cp -r /nfs/mm-isilon/bioinfcore/ActiveProjects/workshop/inputs/* source/inputs 

# On AWS
cp -r ~/ISC_R/inputs/10x_cellranger_filtered_triples source/inputs/
cp -r ~/ISC_R/inputs/prepared_data/phenos.csv source/inputs/prepared_data/
