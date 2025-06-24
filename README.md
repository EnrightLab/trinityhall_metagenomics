## Dana Bao - Summer Research GitHub Repository

Nanopore Chuckling Goat Analysis

The fastq file was obtained directly from the PromethION computer. 
Given the number of files, the fastq files are decompressed, combined and recompressed into 4 barcode files using:

```
cat barcode01/*.fastq.gz > barcode1.fastq.gz
cat barcode02/*.fastq.gz > barcode2.fastq.gz
cat barcode03/*.fastq.gz > barcode3.fastq.gz
cat barcode04/*.fastq.gz > barcode4.fastq.gz
```

Kraken2 was then used to map the nanopore reads against gtdb database.

```
kraken2 --db gtdb --gzip-compressed --threads 20 --output barcode0X.kraken.out --report barcode0X.kraken.report.txt --confidence 0.5 barcodeX.fastq.gz
```
where X represents the number for the combined fastq file.

* `--db gtdb` Use database gtdb, which contains microbial sequences
* `--gzip-compressed` Input reads are in gzip file
* `--threads 20` Use 20 cpus per job
* `--report` Return a summary report
* `--confidence 0.5` Threshold for fraction of k-mers supporting the classification, otherwise considered unclassified

The adapters are removed with super accuracy using raw signal on the PromethION computer.





|table|one|two|
|-----|---|---|
|a    | b | d |

```
somecode here
```


