---
title: "Cell Ranger in Action"
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

## Objective

- By the end of the session you will be able to run `cellranger count` with data from the AGC to transform FASTQ files into a feature-barcode matrix.
- (borrowed heavily from 10x tutorials, but as if with data from AGC)

## Do what 10x does 

### What you get from AGC 

#### Base Data Structure

The basic data package from the AGC includes:

- an **\*.md5** file to validate your data transfer.
- a **DemuxStats_\*.csv** that has some basic metrics about how your samples performed on the sequencer. 
- the **fastq_\*** folder containing the fastq.gz files
- a **README.txt** including details about how your samples were processed.

It looks like the tree below.

```
0000-SR
├── 0000-SR.md5
├── DemuxStats_0000-SR.csv
├── fastqs_0000-SR
│   ├── 0000-SR-1-GEX_S25_R1_001.fastq.gz
│   ├── 0000-SR-1-GEX_S25_R2_001.fastq.gz
│   ├── 0000-SR-2-GEX_S26_R1_001.fastq.gz
│   ├── 0000-SR-2-GEX_S26_R2_001.fastq.gz
│   ├── 0000-SR-3-GEX_S27_R1_001.fastq.gz
│   ├── 0000-SR-3-GEX_S27_R2_001.fastq.gz
│   ├── 0000-SR-4-GEX_S28_R1_001.fastq.gz
│   └── 0000-SR-4-GEX_S28_R2_001.fastq.gz
└── README.txt
```

All you need to run 10x cellranger count are the above 10x sample fastq files and the correct [reference genome file](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads). 
 
If the AGC does basic 10x analysis you also receive a **10x\_analysis\_\*** folder with the following structure and files.

```
0000-SR
|...
└── 10x_analysis_0000-SR
    ├── Sample_0000-SR-1
    │   ├── analysis
    │   ├── cloupe.cloupe
    │   ├── filtered_feature_bc_matrix
    │   ├── filtered_feature_bc_matrix.h5
    │   ├── metrics_summary.csv
    │   ├── molecule_info.h5
    │   ├── possorted_genome_bam.bam
    │   ├── possorted_genome_bam.bam.bai
    │   ├── raw_feature_bc_matrix
    │   ├── raw_feature_bc_matrix.h5
    │   └── web_summary.html
    ├── Sample_0000-SR-2
    │   ├── ...
    ├── Sample_0000-SR-3
    │   ├── ...
    └── Sample_0000-SR-4
        └── ...

```

### How does the AGC generate the 10x_analysis folder and data

***You cannot run the cellranger pipeline easily on a personal computer/laptop because of the [software requirements](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-system-requirements)***, briefly:

> Cell Ranger pipelines run on Linux systems that meet these minimum requirements:
> 
> - 8-core Intel or AMD processor (16 cores recommended). 
> - 64GB RAM (128GB recommended).  
> - 1TB free disk space.  
> - 64-bit CentOS/RedHat 7.0 or Ubuntu 14.04 - See the 10x Genomics OS Support page for details.  
> 
> The pipelines also run on clusters that meet these additional minimum requirements:
> 
> - Shared filesystem (e.g. NFS)
> - Slurm batch scheduling system 

The AGC uses the Great Lake Computing cluster which meets the above requirements.

- **If your lab has signed up for the [RCP](https://arc.umich.edu/umrcp/), you can also get access to the cluster and the available software.** 
- Using the cluster requires an **\*.sbat** file, or slurm batch file. This file allows the slurm scheduler to get the correct resources needed to process your data using cellranger and charge the correct computing account. *If you have questions about this, please ask them in office hours.*
- At the AGC we construct the cellranger command that goes into the *.sbat using one node per sample with 16 cores and 128G of RAM, as recommeded by the software requirements.


#### Constructing the cellranger command (From the 10x Tutorial)
If running interactively on the cluster the first step is to load the cell ranger module with `module load cellranger`.

Print the usage statement to see what is needed to build the command. `cellranger count --help`

The output is similar to the following:

```
cellranger-count
Count gene expression (targeted or whole-transcriptome) and/or feature barcode reads from a single sample and GEM well
 
USAGE:
    cellranger count [FLAGS] [OPTIONS] --id &lt;ID&gt; --transcriptome &lt;PATH&gt;
     
FLAGS:
          --no-bam                  Do not generate a bam file
          --nosecondary             Disable secondary analysis, e.g. clustering. Optional
          --include-introns         Include intronic reads in count
          --no-libraries            Proceed with processing using a --feature-ref but no Feature Barcode libraries
                                    specified with the &#x27;libraries&#x27; flag
          --no-target-umi-filter    Turn off the target UMI filtering subpipeline. Only applies when --target-panel is
                            used
          --dry                     Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop
          --disable-ui              Do not serve the web UI
          --noexit                  Keep web UI running after pipestance completes or fails
          --nopreflight             Skip preflight checks
      -h, --help                    Prints help information
  ...
  
```


To run `cellranger count`, you need to specify an **`--id`**. This can be any string, which is a sequence of alpha-numeric characters, underscores, or dashes and no spaces, that is less than 64 characters. Cell Ranger creates an output directory that is named using this id. This directory is called a &quot;pipeline instance&quot; or pipestance for short.

The `--fastqs` should be a path to the directory containing the FASTQ files. If there is more than one sample in the FASTQ directory, use the **`--sample`** argument to specify which samples to use. This **`--sample`** argument works off of the sample id at the beginning of the FASTQ file name. The last argument needed is the path to the **`--transcriptome`** reference package. Be sure to edit the file paths in the command below.

```
cellranger count --id=0000-SR-1 \
   --fastqs=/nfs/turbo/path/to/0000-SR/0000-SR_fastqs \
   --sample=0000-SR-1\
   --transcriptome=path/to/refereces/cellranger_count/refdata-gex-GRCh38-2020-A`
```

A full-sized dataset can take several hours to complete, but we're using a smaller data set that should process quickly.

The `cellranger count` screen output is similar to the following:

```
/mnt/yard/user.name/yard/apps/cellranger-7.2.0/bin
cellranger count (7.2.0)

cellranger count (7.2.0)
Copyright (c) 2021 10x Genomics, Inc.  All rights reserved.
-------------------------------------------------------------------------------

Martian Runtime - v4.0.6
...
2021-10-15 17:12:42 [perform] Serializing pipestance performance data.
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!
```

When the output of the `cellranger count` command says, “Pipestance completed successfully!”, this means the job is done.
The `cellranger count` pipeline outputs are in the pipestance directory in the outs folder. List the contents of this directory with `ls -1`.

`ls -1 0000-SR-1/outs`

The output is similar to the following:

```
├── analysis 
├── cloupe.cloupe
├── filtered_feature_bc_matrix
├── filtered_feature_bc_matrix.h5
├── metrics_summary.csv
├── molecule_info.h5
├── possorted_genome_bam.bam
├── possorted_genome_bam.bam.bai
├── raw_feature_bc_matrix
├── raw_feature_bc_matrix.h5
└── web_summary.html
```





## QC - key outputs of the quality control report to review 

**QC of small smaple/realistic samples**

***_While cellranger processes we can take a break or discuss what cellranger is doing under the hood if we haven't covered it elsewhere._***

The important outputs we look at to determine the quility of the 10x analaysis include:

* **web_summary.html**
	* **Barcode knee plot**
	* Cell count, UMI count, gene count
	* **Seq saturation**
	* Cell count, UMI count, gene count
* **metrics_summary.csv**



### Web Summary html report. 
The `cellranger count` web summary has Summary and Gene Expression tabs, Similar web summaries are also output from the other cellranger pipelines. 

The run summary from `cellranger count` can be viewed by clicking Summary in the top left tab of the HTML file. The summary metrics describe sequencing quality and various characteristics of the detected cells.

<style data-emotion="css 92repi">.css-92repi{display:block;max-width:100%;}</style><img alt="" class=" css-92repi" src="https://cdn.10xgenomics.com/image/upload/v1646938726/software-support/3p-Single-Cell-GEX/count-web-summary/cr7-count-gexAb-summary1.png"/>
<br/>

#### Summary tab
The estimated number of cells detected, mean reads per cell, and median genes detected per cell are prominently displayed near the top of the page.</p>

Estimated [Sequecing Saturation](https://kb.10xgenomics.com/hc/en-us/articles/115005062366-What-is-sequencing-saturation) is also displayed in the sequencing table.
 
 - The recommended saturation is >60%

Click the `?` icons next to the Sequencing, Mapping, and Cells sections to display information about each metric in the dashboard.</p>

##### The Barcode Rank Plot
The [GEX Barcode Rank Plot](/support/software/cell-ranger/latest/advanced/cr-barcode-rank-plot) under the Cells dashboard shows the distribution of barcode counts and which barcodes were inferred to be associated with cells. 
The y-axis is the number of UMI counts mapped to each barcode and the x-axis is the number of barcodes below that value. 
A steep drop-off is indicative of good separation between the cell-associated barcodes and the barcodes associated with empty partitions. Since barcodes can be associated with cells based on their UMI count or by their RNA profiles, some regions of the graph can contain both cell-associated and background-associated barcodes. 
The color of the graph represents the local density of barcodes that are cell-associated. 

See this [Guided Tour of the Barcode Rank Plot](/support/software/cell-ranger/latest/advanced/cr-barcode-rank-plot) for more details.

#### Gene Expression tab
The automated secondary analysis results can be viewed by clicking the Gene Expression tab in the top left corner. Click the `?` icons next to each section title to display information about the secondary analyses shown in the dashboard.</p>

The t-SNE Projection section shows the data reduced to two dimensions, colored by UMI count (left) or clustering (right). It is a good starting point to explore structure in the data. The projection colored by UMI counts is indicative of the RNA content of the cells and often correlates with cell size - redder points are cells with more RNA in them. For the projection colored by clustering results, select the type of clustering analysis to display from the drop-down button on the upper right (Graph-based by default) - change the category to vary the type of clustering and/or number of clusters (K=2-10) that are assigned to the data.</p>

The Top Features By Cluster table shows which genes are differentially expressed in each cluster relative to all other clusters (Graph-based by default). To find the genes associated with a particular cluster, click the L2FC or p-value column headers associated with a given cluster number to sort the table by a specific cluster</p>
<img alt="" class=" css-92repi" src="https://cdn.10xgenomics.com/image/upload/v1646938726/software-support/3p-Single-Cell-GEX/count-web-summary/cr7-count-gexAb-GEX1.png"/>
<br/>

##### The Sequencing Saturation plot
The Sequencing Saturation plot shows the effect of decreased sequencing depth on sequencing saturation, which is a measure of the fraction of library complexity that was observed. The right-most point on the line is the full sequencing depth obtained in this run.</p>

Similarly, the Median Genes per Cell plot shows the effect of decreased sequencing depth on median genes per cell, which is a way of measuring data yield as a function of depth. The right-most point on the line is the full sequencing depth obtained in this run.</p>
<img alt="" class=" css-92repi" src="https://cdn.10xgenomics.com/image/upload/v1646938726/software-support/3p-Single-Cell-GEX/count-web-summary/cr7-count-gexAb-GEX2.png"/>
<br/>

### Common errors/warnings & remediations
* Saturation & resequencing
* Library prep quality issues

## Outputs used for downstream analysis

[10x output documentation
](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-gex-overview)

### The metrics_summary.csv

Includes much of the same data as the web_summary.html but in csv format.

<table class="css-1pyqnef"><style data-emotion="css eckm6m">.css-eckm6m >tr>*{background-color:#F7F9FC;height:32px;padding:0 16px;}.css-eckm6m >tr>:first-of-type{border-top-left-radius:0.375rem;}.css-eckm6m >tr>:last-of-type{border-top-right-radius:0.375rem;}</style><thead class="css-eckm6m"><tr><style data-emotion="css vmhyu4">.css-vmhyu4{padding:1rem;text-align:left;vertical-align:bottom;}</style><style data-emotion="css yz5ans">.css-yz5ans{-moz-osx-font-smoothing:grayscale;-webkit-font-smoothing:antialiased;box-sizing:border-box;color:#0A2347;font-family:Calibre,sans-serif;font-style:normal;font-weight:500;margin:0;padding:0;font-size:1.125rem;letter-spacing:0.15px;line-height:1.333;padding:1rem;text-align:left;vertical-align:bottom;}</style><th class="css-yz5ans" style="padding:1rem">Column Name</th><th class="css-yz5ans" style="padding:1rem">Description</th></tr></thead><style data-emotion="css 8b62zi">.css-8b62zi >tr>*{background-color:#ffffff;height:32px;padding:0 16px;}.css-8b62zi >:last-of-type>:first-of-type{border-bottom-left-radius:0.375rem;}.css-8b62zi >:last-of-type>:last-of-type{border-bottom-right-radius:0.375rem;}</style><tbody class="css-8b62zi"><tr><style data-emotion="css s9mn2y">.css-s9mn2y{text-align:left;vertical-align:top;}</style><style data-emotion="css 4fk7u9">.css-4fk7u9{-moz-osx-font-smoothing:grayscale;-webkit-font-smoothing:antialiased;box-sizing:border-box;color:#445979;font-family:Calibre,sans-serif;font-style:normal;font-weight:400;margin:0;padding:0;font-size:1.125rem;letter-spacing:0.15px;line-height:1.333;text-align:left;vertical-align:top;}</style><td class="css-4fk7u9" style="padding:1rem">`Estimated Number of Cells`</td><td class="css-4fk7u9" style="padding:1rem">The number of barcodes associated with cell-containing partitions, estimated from the barcode UMI count distribution.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Mean Reads per Cell`</td><td class="css-4fk7u9" style="padding:1rem">The total number of reads divided by the estimated number of cells.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Median Genes per Cell`</td><td class="css-4fk7u9" style="padding:1rem">Median number of read pairs sequenced from the cells assigned to this sample. In case of multiplexing, only cell-associated barcodes assigned exactly one CMO can be assigned to a sample.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Number of Reads`</td><td class="css-4fk7u9" style="padding:1rem">Total number of sequenced reads.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Valid Barcodes`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads with cell-barcodes that match the whitelist.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Sequencing Saturation`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is a ratio where: the denominator is the number of confidently-mapped reads with a valid cell-barcode and valid UMI, and the numerator is the subset of those reads that had a non-unique combination of (cell-barcode, UMI, gene). This metric was called &quot;cDNA PCR Duplication&quot; in versions of Cell Ranger prior to 1.2.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Q30 Bases in Barcode`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of bases with Q-score at least 30 in the cell barcode sequences. This is the i7 index (I1) read for the Single Cell 3&#x27; v1 chemistry and the R1 read for the Single Cell 3&#x27; v2 chemistry.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Q30 Bases in RNA`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of bases with Q-score at least 30 in the RNA read sequences. This is Illumina R1 for the Single Cell 3&#x27; v1 chemistry and Illumina R2 for the Single Cell 3&#x27; v2 chemistry.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Q30 Bases in UMI`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of bases with Q-score at least 30 in the UMI sequences. This is the R2 read for the Single Cell 3&#x27; v1 chemistry and the R1 read for the Single Cell 3&#x27; v2 chemistry.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Reads Mapped to Genome`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads that mapped to the genome.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Reads Mapped Confidently to Genome`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads that mapped uniquely to the genome. If a read mapped to exonic loci from a single gene and also to non-exonic loci, it is considered uniquely mapped to one of the exonic loci.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Reads Mapped Confidently to Intergenic Regions`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads that mapped to the intergenic regions of the genome with a high mapping quality score as reported by the aligner.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Reads Mapped Confidently to Intronic Regions`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads that mapped to the intronic regions of the genome with a high mapping quality score as reported by the aligner.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Reads Mapped Confidently to Exonic Regions`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads that mapped to the exonic regions of the genome with a high mapping quality score as reported by the aligner.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Reads Mapped Confidently to Transcriptome`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads that mapped to a unique gene in the transcriptome with a high mapping quality score as reported by the aligner. The read must be consistent with annotated splice junctions when include-introns=false. These reads are considered for UMI counting.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Reads Confidently Mapped Antisense`</td><td class="css-4fk7u9" style="padding:1rem">Fraction of reads confidently mapped to the transcriptome, but on the opposite strand of their annotated gene. A read is counted as antisense if it has any alignments that are consistent with an exon of a transcript but antisense to it, and has no sense alignments.</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">`Total Genes Detected Median UMI Counts per Cell`</td><td class="css-4fk7u9" style="padding:1rem">The number of genes with at least one UMI count in any cell.</td></tr></tbody></table> 

### Raw\_feature\_bc_matrix / filtered\_feature\_bc\_matrix
	
The `cellranger` pipeline outputs unfiltered (raw) and filtered feature-barcode matrices in two file formats: the Market Exchange Format (MEX), which is described on this page, and Hierarchical Data Format (HDF5), which is described in detail here.

Each element of the feature-barcode matrix is the number of UMIs associated with a feature (row) and a barcode (column):</p>

<style data-emotion="css ayshjd">.css-ayshjd{overflow-x:auto;}</style><div class="css-ayshjd"><style data-emotion="css 1pyqnef">.css-1pyqnef{background-color:#DCE0E6;border-spacing:1px;border-radius:0.375rem;margin-bottom:1.5rem;table-layout:auto;width:100%;}</style><table class="css-1pyqnef"><style data-emotion="css eckm6m">.css-eckm6m >tr>*{background-color:#F7F9FC;height:32px;padding:0 16px;}.css-eckm6m >tr>:first-of-type{border-top-left-radius:0.375rem;}.css-eckm6m >tr>:last-of-type{border-top-right-radius:0.375rem;}</style><thead class="css-eckm6m"><tr><style data-emotion="css vmhyu4">.css-vmhyu4{padding:1rem;text-align:left;vertical-align:bottom;}</style><style data-emotion="css yz5ans">.css-yz5ans{-moz-osx-font-smoothing:grayscale;-webkit-font-smoothing:antialiased;box-sizing:border-box;color:#0A2347;font-family:Calibre,sans-serif;font-style:normal;font-weight:500;margin:0;padding:0;font-size:1.125rem;letter-spacing:0.15px;line-height:1.333;padding:1rem;text-align:left;vertical-align:bottom;}</style><th class="css-yz5ans" style="padding:1rem">Type</th><th class="css-yz5ans" style="padding:1rem">Description</th></tr></thead><style data-emotion="css 8b62zi">.css-8b62zi >tr>*{background-color:#ffffff;height:32px;padding:0 16px;}.css-8b62zi >:last-of-type>:first-of-type{border-bottom-left-radius:0.375rem;}.css-8b62zi >:last-of-type>:last-of-type{border-bottom-right-radius:0.375rem;}</style><tbody class="css-8b62zi"><tr><style data-emotion="css s9mn2y">.css-s9mn2y{text-align:left;vertical-align:top;}</style><style data-emotion="css 4fk7u9">.css-4fk7u9{-moz-osx-font-smoothing:grayscale;-webkit-font-smoothing:antialiased;box-sizing:border-box;color:#445979;font-family:Calibre,sans-serif;font-style:normal;font-weight:400;margin:0;padding:0;font-size:1.125rem;letter-spacing:0.15px;line-height:1.333;text-align:left;vertical-align:top;}</style><td class="css-4fk7u9" style="padding:1rem">Unfiltered feature-barcode matrix</td><td class="css-4fk7u9" style="padding:1rem">Contains every barcode from the fixed list of known-good barcode sequences that has at least one read. This includes background and cell-associated barcodes. <br/>count: `outs/raw_feature_bc_matrix/`<br/> multi: `outs/multi/count/raw_feature_bc_matrix/`</td></tr><tr><td class="css-4fk7u9" style="padding:1rem">Filtered feature-barcode matrix</td><td class="css-4fk7u9" style="padding:1rem">Contains only detected cell-associated barcodes. <br/>count: `outs/filtered_feature_bc_matrix/`<br/> multi: `outs/per_sample_outs/count/sample_filtered_feature_bc_matrix/`</td></tr></tbody></table></div>


For sparse matrices, the matrix is stored in the [Market Exchange Format (MEX)](https://math.nist.gov/MatrixMarket/formats.html). It contains gzipped TSV files with feature and barcode sequences corresponding to row and column indices respectively. For example, the matrices output may look like:</p>

```cd /home/jdoe/runs/sample345/outs
tree filtered_feature_bc_matrix
filtered_feature_bc_matrix
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
0 directories, 3 files
```
Features correspond to row indices. For each feature, the feature ID and name are stored in the first and second column of the (unzipped) `features.tsv.gz` file, respectively. The third column identifies the type of feature, which will be one of `Gene Expression`, `Antibody Capture`, `CRISPR Guide Capture`, `Multiplexing Capture`, or `CUSTOM`, depending on the feature type. Below is a minimal example `features.tsv.gz` file showing data collected for three genes.</p>


```
ENSG00000141510       TP53         Gene Expression
ENSG00000012048       BRCA1        Gene Expression
ENSG00000139687       RB1          Gene Expression
```

Gene Expression</code> data, the ID corresponds to <code class="css-q0dnij">gene_id</code> in the annotation field of the reference GTF. Similarly, the name corresponds to <code class="css-q0dnij">gene_name</code> in the annotation field of the reference GTF. If no <code class="css-q0dnij">gene_name</code> field is present in the reference GTF, gene name is equivalent to gene ID. Similarly, for <code class="css-q0dnij">Antibody Capture</code> and <code class="css-q0dnij">CRISPR Guide Capture</code> data, the <code class="css-q0dnij">id</code> and <code class="css-q0dnij">name</code> are taken from the first two columns of the <a href="/support/software/cell-ranger/latest/analysis/running-pipelines/cr-feature-bc-analysis#feature-ref" class="css-k0f5v edy1dh10">Feature Reference CSV file</a>.</p>

For multi-species experiments, gene IDs and names are prefixed with the genome name to avoid name collisions between genes of different species e.g., GAPDH becomes <code class="css-q0dnij">hg19_GAPDH</code> and Gm15816 becomes <code class="css-q0dnij">mm10_Gm15816</code>.</p>

Barcode sequences correspond to column indices:</p>

```
head filtered_feature_bc_matrices/barcodes.tsv.gz

AAACCCAAGGAGAGTA-1      
AAACGCTTCAGCCCAG-1
AAAGAACAGACGACTG-1
AAAGAACCAATGGCAG-1    
AAAGAACGTCTGCAAT-1    
AAAGGATAGTAGACAT-1    
AAAGGATCACCGGCTA-1    
AAAGGATTCAGCTTGA-1    
AAAGGATTCCGTTTCG-1    
AAAGGGCTCATGCCCT-1
```

Each barcode sequence includes a suffix with a dash separator followed by a number:</p>

<p class="css-d17snu">More details on the barcode sequence format are available in the <a href="/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-bam" class="css-k0f5v edy1dh10">barcoded BAM section</a>.</p>
<p class="css-d17snu">R and Python support the MEX format and sparse matrices can be used for more efficient manipulation.</p>
<p class="css-d17snu">For suggestions on downstream analysis with 3rd party R and Python tools, see the 10x Genomics <a href="/analysis-guides/continuing-your-journey-after-running-cell-ranger" class="css-k0f5v edy1dh10">Analysis Guides</a> resource.</p>
<p class="css-d17snu">The R package Matrix supports loading MEX format data, and can be easily used to load the sparse feature-barcode matrix.</p>

###Converting matrix files to CSV format

Cell Ranger represents the feature-barcode matrix using sparse formats (only the nonzero entries are stored) in order to minimize file size. All of our programs, and many other programs for gene expression analysis, support sparse formats.</p>

However, certain programs (e.g. Excel) only support dense formats (where every row-column entry is explicitly stored, even if it&#x27;s a zero). Here are a few methods for converting feature-barcode matrices to CSV:</p>

###Load matrices into Python
<p class="css-d17snu">The <code class="css-q0dnij">csv</code>, <code class="css-q0dnij">os</code>, <code class="css-q0dnij">gzip</code>, and <code class="css-q0dnij">scipy.io</code> modules can be used to load a feature-barcode matrix into Python.</p>

###mat2csv
<p class="css-d17snu">You can convert a feature-barcode matrix to dense CSV format using the <code class="css-q0dnij">cellranger mat2csv</code> command.</p>

**Important**

**WARNING**: Dense files can be very large and may cause Excel to crash, or even fail in mat2csv if your computer doesn&#x27;t have enough memory. For example, a feature-barcode matrix from a human reference (~33k genes) with ~3k barcodes uses at least 200MB of disk space. Our 1.3 million mouse neuron dataset, if converted to this format, would use more than 60GB of disk space. Thus, while you can use mat2csv for small datasets, we strongly recommend using R or Python (as shown in the sections above) to examine these matrix files.

This command takes two arguments - an input matrix generated by Cell Ranger (either an HDF5 file or a MEX directory), and an output path for the dense CSV. For example, to convert a matrix from a pipestance named <code class="css-q0dnij">sample123</code> in the current directory, either of the following commands would work:</p>

```
# Convert from MEX
cellranger mat2csv sample123/outs/filtered_feature_bc_matrix sample123.csv

# Or, convert from HDF5
cellranger mat2csv sample123/outs/filtered_feature_bc_matrix.h5 sample123.csv
```
* H5 vs files

[intro to HD5 format](https://portal.hdfgroup.org/documentation/index.html)

[10x HDF docs](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-molecule-info)

The cellranger pipeline outputs an HDF5 file containing per-molecule information for all molecules that contain a valid barcode, a valid UMI, and were assigned with high confidence to a gene or Feature Barcode. This HDF5 file contains data corresponding to the observed molecules, as well as data about the libraries and feature set(s) used (general information about the HDF5 file format available here). This file is called molecule_info.h5 in cellranger count and sample_molecule_info.h5 in cellranger multi outputs.

```
(root)
    ├─ barcode_idx
    ├─ barcode_info	[HDF5 group]
    │   ├─ genomes
    │   └─ pass_filter
    ├─ barcodes
    ├─ count
    ├─ feature_idx
    ├─ features	[HDF5 group]
    │   ├─ _all_tag_keys
    │   ├─ target_sets [for Targeted Gene Expression or Fixed RNA Profiling]
    │   │    └─ [target set name]
    │   ├─ feature_type
    │   ├─ genome
    │   ├─ id
    │   ├─ name
    │   ├─ pattern [Feature Barcode only]
    │   ├─ read [Feature Barcode only]
    │   └─ sequence [Feature Barcode only]
    ├─ gem_group
    ├─ library_idx
    ├─ library_info
    ├─ metrics_json [HDF5 dataset; see below]
    ├── probe_idx          ---------------------|
    ├── probes [HDF5 group]                     |
    │   ├── feature_id                          | [For Fixed RNA Profiling, Cell Ranger v7.1+]
    │   ├── feature_name                        |
    │   ├── probe_id                            |
    │   └── region         ---------------------|
    ├─ umi
    └─ umi_type
```

You can examine the contents of the H5 file using software such as HDFView or the h5dump command, as demonstrated below to show the file contents of the entire H5 object:

```h5dump -n molecule_info.h5

    HDF5 "molecule_info.h5" {
    FILE_CONTENTS {
    group      /
    dataset    /barcode_idx
    group      /barcode_info
    dataset    /barcode_info/genomes
    dataset    /barcode_info/pass_filter
    dataset    /barcodes
    dataset    /count
    dataset    /feature_idx
    group      /features
    dataset    /features/_all_tag_keys
    dataset    /features/feature_type
    dataset    /features/genome
    dataset    /features/id
    dataset    /features/name
    dataset    /gem_group
    dataset    /library_idx
    dataset    /library_info
    dataset    /metrics_json
    dataset    /umi
    dataset    /umi_type
    }
    }
```

More details can be found at the following links:

[intro to HD5 format](https://portal.hdfgroup.org/documentation/index.html)

[10x HDF docs](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-molecule-info)

<br/>
<br/>
<hr/>
| [Previous lesson](00A-OrientingOnScRNASeq.html) | [Top of this lesson](#top) | [Next lesson](01-GettingStarted.html) |
