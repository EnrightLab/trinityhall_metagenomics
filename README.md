## Dana Bao - Summer Research GitHub Repository

**Nanopore Chuckling Goat Analysis**

The fastq file was obtained directly from the PromethION computer. 
Given the number of files, the fastq files were decompressed, combined and recompressed into 4 barcode files:

```
cat barcode0X/*.fastq.gz > barcodeX.fastq.gz
```
where X represents the number for the combined fastq file.

Kraken2 was then used to map the nanopore reads against gtdb and core_nt database.

```
kraken2 --db gtdb --gzip-compressed --threads 20 --output barcode0X.kraken.out --report barcode0X.kraken.report.txt --confidence 0.5 barcodeX.fastq.gz
```
```
kraken2 --db kraken_db --gzip-compressed --threads 20 --output barcode0X.kraken.out --report barcode0X.kraken.report.txt --confidence 0.5 barcodeX.fastq.gz
```

* `--db gtdb` or `--db kraken_db`: Use database gtdb/core_nt, which contains microbial sequences/full NCBI's database sequences
* `--gzip-compressed`: Input reads are in gzip file
* `--threads 20`: Use 20 cpus per job
* `--report`: Return a summary report
* `--confidence 0.5`: Threshold for fraction of k-mers supporting the classification, otherwise considered unclassified

Sample Kraken outputs are as follow:
```
71.37  16512460        16512460        U       0       unclassified
 28.63  6624973 17337   R       1       root
 28.56  6607634 266707  D       5839      d__Bacteria
 27.37  6332637 347     P       28625       p__Firmicutes
 27.37  6332290 55784   C       28626         c__Bacilli
 27.12  6274870 63525   O       39128           o__Lactobacillales
 26.83  6208064 531965  F       107261            f__Lactobacillaceae
 20.89  4833090 1588432 G       122880              g__Lactobacillus
 12.02  2781778 0       S       153258                s__Lactobacillus kefiranofaciens
 12.02  2781778 2781778 S1      153265                  RS_GCF_900103655.1
  1.97  456293  0       S       128401                s__Lactobacillus helveticus
  1.97  456293  456293  S1      128411                  RS_GCF_000160855.1
  0.01  2728    0       S       177813                s__Lactobacillus ultunensis
  0.01  2728    2728    S1      177815                  RS_GCF_001436305.1
```

The Kraken report was used to perform abundance estimation using Bracken: 
```
est_abundance.py -i barcode0X.kraken.report.txt -k kraken_db/database75mers.kmer_distrib -l S -t 10 -o barcode0X.bracken.txt
```
* `--i`: Input kraken report file name
* `--k`: The k-mer distribution file used for the previous kraken analysis
* `--l S`: Species level abundance estimation
* `--t 10`: Use 10 cpus per job
* `--o`: Output bracken file name

Sample Bracken outputs are as follow:
```
name    taxonomy_id     taxonomy_lvl    kraken_assigned_reads   added_reads     new_est_reads   fraction_total_reads
s__Lactobacillus kefiranofaciens        153258  S       2781778 0       2781778 0.69819
s__Lactobacillus helveticus     128401  S       456293  0       456293  0.11452
s__Lactobacillus ultunensis     177813  S       2728    0       2728    0.00068
s__Lactobacillus gallinarum     149202  S       2043    0       2043    0.00051
s__Lactobacillus crispatus      122881  S       1039    0       1039    0.00026
s__Lactobacillus amylovorus     141612  S       630     2       632     0.00016
s__Lactobacillus amylolyticus   156383  S       55      0       55      0.00001
s__Lactobacillus kitasatonis    186364  S       54      0       54      0.00001
s__Lactobacillus acetotolerans  139478  S       35      0       35      0.00001
s__Lentilactobacillus kefiri    157081  S       694001  17134   711135  0.17849
s__Lentilactobacillus parakefiri        166863  S       7826    0       7826    0.00196
```
For abundance visualization, Krona was used. However, the Krona in-built taxonomy uses NCBI taxonomy id instead of gtdb taxnomy id, so the conversion between the two was required for bracken outputs generated against gtdb database. Bracken outputs derived from searching against core_nt uses NCBI taxonomy id therefore do not need conversion. Using version r207 of gtdb release as an example, first download the bacterial and archaean metadata:

```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/ar53_metadata_r207.tar.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/bac120_metadata_r207.tar.gz
```
Then construct a full metadata file using `bac120` and `ar53` via python. As the metadata lack gtdb taxonomy id, species name were isolated instead to prepare for later mapping. Species names were extract from `gtdb_taxonomy`column:

```
import pandas as pd

bac_metadata = pd.read_csv("$PATH/bac120_metadata_r207.tsv", sep="\t")
ar_metadata = pd.read_csv("$PATH/ar53_metadata_r207.tsv", sep="\t")

metadata = pd.concat([bac_metadata, ar_metadata], ignore_index=True)
metadata["s_name"]= (metadata["gtdb_taxonomy"].str.extract(r'(s__.+)$').fillna("").squeeze())

metadata.to_csv("$PATH/all_microbes_metadata.tsv", sep="\t", index=False)
```

The full metadata was then used to map gtdb taxonomy id from Bracken output onto NCBI taxonomy id by corresponding species names via python:
```
import pandas as pd

metadata = pd.read_csv("$PATH/all_microbes_metadata.tsv", sep="\t", dtype={"ncbi_taxid":"Int64"})
inputfile = input("Enter input bracken file full path:")
input_df = pd.read_csv(inputfile, sep="\t")

species_to_taxid = dict(zip(metadata["s_name"], metadata["ncbi_taxid"]))
input_df['ncbi_taxid'] = input_df['name'].map(species_to_taxid)

input_df["taxonomy_id"] = input_df["ncbi_taxid"]
input_df = input_df.drop(columns=["ncbi_taxid"])
input_df["taxonomy_id"] = input_df["taxonomy_id"].fillna(0).astype(int)

outputfile = inputfile.replace(".bracken.txt", ".bracken.ncbi.txt")
input_df.to_csv(outputfile, sep="\t", index=False)
```

The output `barcode0X.bracken.ncbi.txt`, along with `barcode0X.bracken.core.txt`, were used as input for generating Krona plot:
```
ktImportTaxonomy -t 2 -m 6 -o barcode0X.html barcode0X.bracken.ncbi.txt
```
```
ktImportTaxonomy -t 2 -m 6 -o barcode0X.html barcode0X.bracken.core.txt
```

* `--t 2`: Use the 2nd column of Bracken output for taxonomy id search, which contains NCBI taxonomy id of the species
* `--m 6`: Use the 6th column of Bracken output for abundance score, which contains `new_est_reads`, which is the bayesian restimation of abundance based on Kraken output

Sample Krona Plot is as follow:
<img width="892" alt="barcode01" src="https://github.com/user-attachments/assets/a16007aa-3229-4ffa-8671-829f3c6b192c" />

Analysis from Krona plot identified 5 major species of microbial organisms in the samples. They were Lactobacillus kefiranofaciens, Lactobacillus helveticus, Lentilactobacillus kefiri, Kluyeveromyces marxianus and Pichia kudriavzevii, respectively. 

In order to further analyze these species, their specific genome sequences needs to be created using nanopore reads. To filter out reads of interest, a python script `find_species_hits.py` was used to identify read ids corresponding to the species of interest from kraken output file, then filtering out the corresponding sequence reads from the fastq file. This was done for all 4 samples, and the resulting sequences for each species were combined into a fasta file.

```
python find_species_hits.py barcode0X.kraken.core.out
```

Reference genomes were downloaded from Genbank: 
|Species                                  |Genbank Accession|
|-----------------------------------------|-----------------|
|Lactobacillus kefiranofaciens strain 1207| NZ_CP061341     |
|Lactobacillus helveticus strain TCI357   | NZ_CP151471     |
|Lentilactobacillus kefiri strain DH5     | NZ_CP029971     |
|Kluyveromyces marxianus DMKU3-1042       | NC_036025       |
|Pichia kudriavzevii                      | NC_042506       |

Alignment of metagenomic data to reference genomes is done through blast. To produce desired output format, `blaster.py` was used:

```
python blaster.py
```
The code requires file `targets.txt` containing species name and corresponding database listed in tsv format to run. 

Sample output is as follow: 

```
6b434867-92f3-46f4-91c3-60e25aa370fa    NZ_CP061341.1   100.000 70.85% [175 / 247]      2146803 2146629 71      245
#
#[aln]  CGGCTCAAGAAGGTAAAATCTTGCGAAATGGCTTGGCTACTGCAATCGTAGGTCGACCAAATGTGGGAAAATCATCCTTG
#[aln]  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#[aln]  CGGCTCAAGAAGGTAAAATCTTGCGAAATGGCTTGGCTACTGCAATCGTAGGTCGACCAAATGTGGGAAAATCATCCTTG
#
#[aln]  CTTAATTATTTAACGCAAAGTGATAAAGCAATTGTAACCGATGTCGCTGGAACTACGCGTGACACCTTGGAAGAATATGT
#[aln]  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#[aln]  CTTAATTATTTAACGCAAAGTGATAAAGCAATTGTAACCGATGTCGCTGGAACTACGCGTGACACCTTGGAAGAATATGT
#
#[aln]  ATCTGTAAAAGGCGT
#[aln]  |||||||||||||||
#[aln]  ATCTGTAAAAGGCGT
#
```

To visualize the comparison between the sequenced strain and reference genome, output from Blast was used to generate Circos plots. Files that are necessary for Circos was created based on Blast output using `make_karyotype_gb.py`. GC skew data was generated using `gcskew.py` (Jennifer Lu, jlu26@jhmi.edu). 

```
python gcskew.py -i species_name_reference_genome.fna -o species_name_skew.txt
```
```
make_karyotype_gb.py species_name.gbff species_name.metagenomic.fasta.blasthits.txt species_name.skew.txt
```

As both `label_forward.txt` and `label_reverse.txt` contain annotations with space and would thus cause Circos script to fail, spaces were substituted with underscores using:

```
sed -i 's/ /_/g' labels_forward.txt
```
```
sed -i 's/ /_/g' labels_reverse.txt
```

The files were then passed to Circos using `circos.conf` as the configuration reference. 
```
circos -conf circos.conf
```

Sample Circos plots are as follow:
![Lactobacillus_kefiranofaciens_circos](https://github.com/user-attachments/assets/d6f9cd17-5996-4fd0-8585-336a4b06c1c1)

For higher processing speed and visualization via Integrative Genomics Viewer, minimap2 was also used to align the read to the reference genome, which is more efficient when applied to long Oxford Nanopore reads. 

```
minimap2 -t 10 -a -x map-ont species_name_reference_genome.fna species_name.metagenomic.fasta > species_name.minimap2.sam
```
* `--a`: Output as SAM format, default is PAF format
* `--t 10`: Use 10 cpus per job
* `--x map-ont`: Specifies input as Oxford Nanopore reads

For visualization using IGV, binary bam files instead of sam files were necessary. This is created using Samtools:

```
samtools view -Sb species_name.minimap2.sam > species_name.bam
```

This produced unsorted bam file without indexing. To create sorted, indexed bam file for IGV, samtools were again used:

```
samtools sort -o species_name.sorted.bam species_name.bam
```
```
samtools index species_name.sorted.bam
```
This creates `species_name.sorted.bam.bai` index file, which allow IGV to locate the alignments in specific regions. `species_name_reference_genome.fna`, `species_name.sorted.bam` and `species_name.sorted.bam.bai` are all loaded into IGV, with sample visualization region as follow:
<img width="1367" height="808" alt="Screenshot 2025-07-29 at 11 20 30" src="https://github.com/user-attachments/assets/a7710b32-379d-46fa-b5b7-38318d4f2cc8" />
The adapters were removed with super accuracy using raw signal on the PromethION computer. The new sequence files acquired were again recompressed into `bracodeX_super_trimmed.fastq.gz`, then handled similarly with Kraken and Bracken.

In addition, porechop was used to trim the adapters. To prevent reaching memory limit, the porechop was run on each individual fastq files instead of on compressed fastq files using perl:

```
#!/usr/bin/perl

foreach $file (@ARGV){


        $output_file=$file;
        $output_file=~ s/.fastq.gz/.chopped.fastq.gz/g;

        system("/data/user_scripts/Porechop/porechop-runner.py -i $file -o $output_file");
}
porechop.pl (END)
```

**Nanopore Soil Analysis**

3 soil samples acquired in Trinity Hall, Cambridge and 1 supersoil sample was analyzed in this experiment. 

Similar to what is done with Chuckling Goat, the data were first processed via Kraken2 to acquire mapping against gtdb database. The confidence level was tested between 0.1 and 0.5 at an incremented step of 0.1, and while higher confidence level yield high amount of unclassified reads, downstream Bracken analysis and Krona plots produces similar result between 0.2 and 0.5 confidence level. Confidence level test results could be found in `barcode13_confidence.microb.html`, result for all 4 soil samples analyzed at confidence level 0.2 could be found in `barcode_confidence02.microb.html`. 

Genus level abundance was also analyzed in this study using Bracken:
```
est_abundance.py -i barcode13.kraken.microb.report.txt -k /mnt/cgs-fs7.hmg.path.private.cam.ac.uk-cgs-fs7/Anton/hazal/kraken_db/database75mers.kmer_distrib -l G -t 10 -o barcode13.G.bracken.microb.txt
```
* `--l G`: Genus level abundance estimation

To visualize the abundance level for both species and genus level, the data acquired at 0.1 confidence level was plotted as stacked barchart for all 4 samples using python script `abundance_plot.ipynb`. Output plots are as follow:
<img width="6291" height="2953" alt="top20_species_confidence02" src="https://github.com/user-attachments/assets/0d72864e-aec0-46da-bdc8-0cb59093ec41" />
<img width="6129" height="2953" alt="top20_genus_confidence02" src="https://github.com/user-attachments/assets/c7a0efc1-6b85-4d45-ab3a-6eea4f5d3e61" />

To quantify the species diversity of each sample and compare their relative abundance, Chao1 index, Shannon index and Simpson index were calculated from the Kraken2 report using direct read numbers for each species (S) or strain/subspecies (S1). 
The Chao1 index is calculated as:  
$S_{chao1} = S_{obs} + \frac{F_1^2}{2F_2}$  
where $S_{obs}$ is the observed species number, $F_1$ is the number of species observed only once, and $F_2$ is the number of species observed twice.

The Shannon index is calculated as:  
$H' = -\sum_{i=1}^{S} p_i \log p_i$  
where $p_i$ is the percentage of individuals belonging to species $i$.

The Simpson index is calculated as:  
$D = \sum_{i=1}^{S} p_i^2$  
where $p_i$ is the percentage of individuals belonging to species $i$.

Script used for calculation and plotting of indexes may be found in diversity_index_calculation.ipynb. Output plots are as follow:

<img width="420" height="293" alt="chao1_confidence02" src="https://github.com/user-attachments/assets/099418c2-2aac-40e4-a8e4-5b9ef6078771" />
<img width="410" height="293" alt="shannon_confidence02" src="https://github.com/user-attachments/assets/96b0d849-54fe-48ef-bcd6-c524a407a145" />
<img width="424" height="293" alt="simpson_confidence02" src="https://github.com/user-attachments/assets/e575fb32-57b2-418f-9200-61a12f0d3932" />

Chao1 index indicates species richness as reflected by number of species and gives more weight to rare species. It can be interpreted from the plot that Sample 1 contains highest number of species, while Sample 4 contains the lowest. 

Shannon index indicates species eveness. It can be interpreted from the plot that Sample 4 has the highest Shannon index, thus having the most even distribution between species. Sample 4 has the lowest index and therefore could have a skewed distribution. 

Simpson index indicates dominance. The results agrees with Shannon index results, with Sample 4 having a low index and thus a low dominance, while Sample 2 having a high index and thus a high dominance of a few species. 
