# Dana Bao - Summer Research GitHub Repository

## Overview

This repository contains bioinformatics analyses of Oxford Nanopore metagenomic sequencing data from Chuckling Goat kefir and soil samples collected in Cambridge. The project characterizes microbial community composition and performs reference-guided strain comparison for dominant taxa, including long-read mapping visualisation and draft genome assembly from species-filtered reads. 
It also tackles practical challenges originating from Oxford Nanopore high error rate, cross-database taxonomic incompatibility, and reference bias in alignments.
The project aims to provide reference and guidance for future projects that involve similar analysis in the group. 

## Data and Experimental Context

Raw FASTQ files were generated directly from a PromethION sequencer and demultiplexed by barcode through Dorado-based calling. 
Dorado demultiplexing errors exists, but are expected to have limited impact on kefir samples, as the samples are duplicates and share the same underlying community structure. 
Contrastingly, demultiplexing accuracy is more consequential for the soil samples, which represent 4 biologically distinct communities. 
This factor is therefore considered in downstream analysis. 
Due to the large number of output files, reads were concatenated at the barcode level prior to downstream analysis.

Two datasets are analysed in this repository:
* Kefir samples dominated by lactic acid bacteria and yeast
* Soil samples with higher taxonomic diversity and lower dominance

## Key Bioinformatics Challenges

The analysis in this project were shaped by several bioinformatics challenges:
* High base-calling error rates of Oxford Nanopore sequencing, affecting k-mer–based taxonomic classification
* Trade-offs between sensitivity and specificity when setting taxonomic classification confidence thresholds
* Incompatibility between GTDB and NCBI taxonomic identifiers across commonly used tools
* Reference bias when extracting species-specific reads from metagenomic data
* Limitations of reference genomes for strain-level comparison and assembly

These challenges therefore affects both tool selection and parameter optimization throughout the project.

## Computational Strategy and Decisions

### Taxonomic classification and confidence thresholding
Kraken2 was used for taxonomic classification against both GTDB and NCBI core_nt databases. For Kraken2 to consider making a taxonomic call, the proportion of k-mers in the read that support the chosen taxon among all informative k-mers in the read must exceed a chosen threshold, namely the confidence threshold. To assess the robustness of classification results, confidence thresholds ranging from 0.1 to 0.5 were compared and evaluated.
While higher confidence thresholds significantly increased the proportion of unclassified reads, downstream abundance estimation using Bracken and visualisation via Krona produced largely consistent community profiles across confidence values ≥0.2. This suggested that community-level patterns were robust to moderate changes in classification stringency, despite read-level uncertainty. 
However, a confidence threshold of 0.1 significantly altered the final community-level profile for soil samples, revealing the risk of false-positive classifications and overestimation of diversity when classification stringency is too low, particularly in situations where community complexity is high and many low-abundance taxa are present. Correspondingly, this effect was less pronounced in the kefir samples, which are dominated by a small number of taxa, making the community profile more robust to low-confidence misclassifications.

### Database choice
Results generated from GTDB-based and NCBI-based workflows were compared. Although NCBI-based databases provide broad coverage across the tree of life, this project focused on microbial community profiling. GTDB offers a genome-based taxonomy curated specifically for bacteria and archaea, which keeps taxonomic outputs aligned with microbiome analysis and reduces distraction from macroorganism signals observed in broad NCBI classifications (e.g., mammalian assignments in soil/kefir datasets). 
In addition, GTDB-based classification produced more informative microbial labels for downstream community profiling, whereas NCBI-based classification more frequently returned low-resolution annotations (e.g., “uncultured bacterium”) for soil sample reads. 

### Cross-database taxonomic harmonisation
A major practical challenge was the incompatibility between local GTDB taxonomy, which was used for better microbiome classification, and NCBI taxonomy, which was required by Krona visualization. To address this, GTDB metadata were parsed and a custom mapping between GTDB species names and NCBI taxonomic identifiers was constructed in Python.
This approach necessarily relies on species-name matching, which may introduce ambiguity in cases of inconsistent nomenclature or unresolved taxonomy. However, it enabled cross-tool compatibility and downstream visualization while maintaining awareness of its limitations. 
In this study, only a small number of low abundance taxa lacked an explicit GTDB to NCBI mapping. These entries were excluded from Krona inputs and had little affect on the community-level patterns.

### Species-level read extraction from metagenomic data
To enable strain-level analysis for kefir sample, reads classified to dominant species were extracted from the metagenomic data based on Kraken output. This step introduces potential misassignment due to classification errors and shared k-mers between closely related taxa, particularly with the higher error rate of Nanopore data.
Nevertheless, this approach provide more specificity than direct metagenome alignment, and the issues arise were partially addressed with downstream alignment, assembly, and visualization choices.

### Reference-guided alignment and visualisation
Extracted reads were aligned to reference genomes using both BLAST and minimap2. 
Minimap2 was favored for large-scale, long-read alignment, and can better accommodate high error rate Nanopore data. 
In this study, BLAST alignments required ~1–3 days, whereas minimap2 completed in ~30 minutes.  
Alignments were visualized using IGV and Circos to assess genome coverage patterns and potential structural divergence, which help reveal misassignments and reference biases. 
No signal representing significant misassignment was observed in this study, suggesting that species-specific read extraction is relatively suitable for dominant species analysis. 

### Strain-level de novo assembly and polishing 
Miniasm was selected for de novo assembly due to its suitability for Nanopore long reads and computational efficiency. As Miniasm constructs assembly graphs without consensus polishing, Racon was used to iteratively polish assemblies using overlap alignments generated by minimap2.
Assembly parameters were adjusted to balance computational feasibility against alignment sensitivity, particularly given the large number of reads and high memory demands.

## Ecological Analysis
Community diversity of soil samples was quantified using Chao1, Shannon, and Simpson indices calculated directly from classified read counts at species and strain levels. These metrics provide complementary perspectives on richness, evenness, and dominance across samples.
This approach assumes proportionality between read counts and organism abundance. This assumption is undermined by genome size differences and sequencing bias, but remains informative for relative comparisons between samples.

## Limitations
*Demultiplexing errors may lead to misassigned reads between barcodes
*Taxonomic classification accuracy may be affected by Nanopore error rates and reference database completeness
*Species-name–based taxonomy mapping between database may introduce ambiguity
*Abundance estimates based on read counts may not directly reflect true abundance
*Strain-level assemblies may be fragmented and require further polishing and validation, ideally supported by short-read data

## Detailed Methods
For a full step-by-step description of the computational workflow, see `pipeline_overview.md`. 
