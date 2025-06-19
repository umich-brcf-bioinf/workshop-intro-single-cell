## Setup input directory, assuming you're at repo root
mkdir -p source/inputs/prepared_data

## Pulls lastest version of inputs (with updated CellRanger) - 1.28.25
#cp -r /nfs/mm-isilon/bioinfcore/ActiveProjects/workshop/inputs/* source/inputs 

## On AWS - file versions depend on syncing source files from core storage if any changes since previous workshop
cp -r ~/ISC_R/inputs/10x_cellranger_filtered_triples source/inputs/
cp -r ~/ISC_R/inputs/prepared_data/*.csv source/inputs/prepared_data/
cp -r ~/ISC_R/inputs/prepared_data/rdata source/inputs/prepared_data/
