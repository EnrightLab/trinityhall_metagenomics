# Dana Bao - Summer Research GitHub Repository

## Overview

This repository is a bioinformatics project focusing on metagenomic analysis of Oxford Nanopore sequencing data generated from Chuckling Goat kefir and soil samples collected in Cambridge. The project involve analysis on microbial community composition, reference-guided strain comparison, long-read alignment visualisation and genome assembly.  
It also addresses some challenges originating from Oxford Nanopore’s high error rate, cross-database taxonomic incompatibility, and reference bias in alignments.  
The project aims to provide reference and guidance for future projects that involve similar analysis in the group.  

## Data and Experimental Context

FASTQ files were obtained directly from the PromethION computer and demultiplexed through Dorado-based calling.  
Dorado demultiplexing errors exists, but are expected to have limited impact on kefir samples, as the samples are duplicates and share the same underlying community structure.  
Contrastingly, this may have more impact on the soil samples, which originate from four different communities.  
This factor is therefore considered in downstream analysis.  

Two datasets are analysed: 
* Kefir samples dominated by lactic acid bacteria and yeasts
* Soil samples with higher taxonomic diversity and lower dominance

## Key Bioinformatics Challenges

* High base-calling error rates of Oxford Nanopore sequencing affects k-mer–based taxonomic classification
* Balancing sensitivity and specificity with taxonomic classification confidence thresholds
* Incompatibility between GTDB and NCBI taxonomic identifier
* Reference bias when extracting species-specific reads from metagenomic data
* Limitations of using reference genomes for strain-level comparisons and assemblies

These challenges therefore affects both tool selection and parameter optimization throughout the project. 

## Computational Strategy and Decisions

### Taxonomic classification and confidence threshold
Kraken2 was used for taxonomic classification against both GTDB and NCBI core_nt databases. For Kraken2 to consider making a taxonomic call, the proportion of k-mers in the read that support the taxon among all informative k-mers must exceed a chosen threshold, namely the confidence threshold. To assess its impact on results, confidence thresholds ranging from 0.1 to 0.5 were compared.  
While higher confidence thresholds significantly increased the proportion of unclassified reads, downstream abundance estimation using Bracken and visualisation via Krona produced largely consistent community profiles for confidence values from 0.2 to 0.5. This suggested that community-level patterns were mostly unaffected by moderate changes in threshold stringency, despite read-level differences.  
However, a confidence threshold of 0.1 significantly altered the final community-level profile for soil samples, revealing the risk of false-positive and diversity overestimation when classification stringency is too low, particularly in situations where community complexity is high and many low-abundance taxa are present. Correspondingly, the kefir samples, which are dominated by a small number of taxa, still remain largely unaffected.  

### Database choice
Results generated from GTDB-based and NCBI-based workflows were compared. As GTDB targets specifically bacteria and archaea, it keeps taxonomic outputs on microbiome and reduces distraction from macroorganism signals observed in broad NCBI classifications (e.g., mammalians in soil/kefir samples).  
In addition, GTDB-based classification assigns more specific taxonomic labels for downstream community profiling, whereas NCBI-based classification more frequently returned low-resolution labels (e.g., “uncultured bacterium”) for soil sample reads.  

### Cross-database taxonomic harmonisation
There is an incompatibility between local GTDB taxonomy, which was used for better microbiome classification, and local NCBI taxonomy, which was used by Krona visualization. To address this, a script that mapped GTDB taxaID back to NCBI taxaID using species name was constructed.  
This approach relies on species-name, which may introduce ambiguity in cases of inconsistent nomenclature or unresolved taxonomy. However, it enabled cross-tool compatibility and downstream visualization.  
In this study, only a small number of low abundance taxa lacked an explicit GTDB to NCBI mapping. They thus had little effect on the community-level patterns.  

### Species-level read extraction from metagenomic data
To enable strain-level analysis for kefir sample, reads classified to dominant species were extracted from the metagenomic data based on Kraken output. This introduces potential misassignment due to classification errors and shared k-mers between closely related taxa, particularly with the higher error rate of Nanopore data.  
Nevertheless, this approach provides higher specificity than direct metagenomic alignment, and the issues were partially addressed with downstream alignment, assembly, and visualization choices.  

### Reference-guided alignment and visualisation
Extracted reads were aligned to reference genomes using both BLAST and Minimap2.  
Minimap2 can better accommodate high error rate Nanopore data, and has a higher alignment efficiency. In this study, BLAST alignments required ~1–3 days, while minimap2 required only ~30 minutes.  
Alignments were visualized using IGV and Circos. This helped reveal misassignments and reference biases.  
There was no signal representing significant misassignment in this study, suggesting that species-specific read extraction is relatively suitable for dominant species analysis with sufficient sequencing depth.  

### Strain-level de novo assembly and polishing
Miniasm was selected for de novo assembly for Nanopore long reads suitability and high computational efficiency. As Miniasm does not contain polishing, Racon was used to iteratively polish assemblies using overlap alignments generated by minimap2.  
Assembly parameters were adjusted to balance computational feasibility and alignment sensitivity, which is essential given the large number of reads and high memory demands.  

## Ecological Analysis
Community diversity of soil samples was quantified using Chao1, Shannon, and Simpson indices. These metrics provide complementary perspectives on richness, evenness, and dominance across samples.  
This approach assumes that read counts are proportional to organism abundance. This is undermined by genome size differences and sequencing bias, but remains informative for relative comparisons between samples.  

## Limitations
*Demultiplexing errors may lead to misassigned reads between barcodes
*Taxonomic classification accuracy may be affected by Nanopore error rates and reference database completeness
*Taxonomy mapping based on species name may introduce ambiguity
*Abundance estimates based on read counts may not directly reflect true abundance
*Strain-level assemblies may be fragmented and require further polishing and validation, ideally supported by short-read data

## Pipeline Overview
For a full description of the computational workflow, see `pipeline_overview.md`.  