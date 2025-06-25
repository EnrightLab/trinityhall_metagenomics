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
where X represents the number for the combined fastq file.

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
where X represents the number for the combined fastq file.

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



