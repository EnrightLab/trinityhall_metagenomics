## Dana Bao - Summer Research GitHub Repository

**Nanopore Chuckling Goat Analysis**

The fastq file was obtained directly from the PromethION computer. 
Given the number of files, the fastq files are decompressed, combined and recompressed into 4 barcode files:

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
For abundance visualization, Krona is used. However, the Krona in-built taxonomy follows NCBI taxonomy id instead of gtdb taxnomy id, so the conversion between the two is required. Using version r207 of gtdb release as an example, first download the bacterial and archaean metadata:

```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/ar53_metadata_r207.tar.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/bac120_metadata_r207.tar.gz
```
Then construct a full metadata file using `bac120` and `ar53` via python. As the metadata lack gtdb taxonomy id, species name are isolated instead to prepare for later mapping. Species names are extract from `gtdb_taxonomy`column:

```
import pandas as pd

bac_metadata = pd.read_csv("$PATH/bac120_metadata_r207.tsv", sep="\t")
ar_metadata = pd.read_csv("$PATH/ar53_metadata_r207.tsv", sep="\t")

metadata = pd.concat([bac_metadata, ar_metadata], ignore_index=True)
metadata["s_name"]= (metadata["gtdb_taxonomy"].str.extract(r'(s__.+)$').fillna("").squeeze())

metadata.to_csv("$PATH/all_microbes_metadata.tsv", sep="\t", index=False)
```

The full metadata is then used to map gtdb taxonomy id from Bracken output onto NCBI taxonomy id by corresponding species names via python:
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

The output `barcode0X.bracken.ncbi.txt` are used as input for generating Krona plot:
```
ktImportTaxonomy -t 2 -m 6 -o barcode0X.html barcode0X.bracken.ncbi.txt
```

* `--t 2`: Use the 2nd column of Bracken output for taxonomy id search, which contains NCBI taxonomy id of the species
* `--m 6`: Use the 6th column of Bracken output for abundance score, which contains `new_est_reads`, which is the bayesian restimation of abundance based on Kraken output

Sample Krona Plot is as follow:
<img width="1189" alt="barcode01" src="https://github.com/user-attachments/assets/a16007aa-3229-4ffa-8671-829f3c6b192c" />

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

|table|one|two|
|-----|---|---|
|a    | b | d |



