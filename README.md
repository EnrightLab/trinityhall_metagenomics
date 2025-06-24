## Dana Bao - Summer Research GitHub Repository

**Nanopore Chuckling Goat Analysis**

The fastq file was obtained directly from the PromethION computer. 
Given the number of files, the fastq files are decompressed, combined and recompressed into 4 barcode files:

```
cat barcode01/*.fastq.gz > barcode1.fastq.gz
cat barcode02/*.fastq.gz > barcode2.fastq.gz
cat barcode03/*.fastq.gz > barcode3.fastq.gz
cat barcode04/*.fastq.gz > barcode4.fastq.gz
```

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

The adapters were removed with super accuracy using raw signal on the PromethION computer.

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



