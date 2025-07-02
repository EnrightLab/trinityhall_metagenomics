import sys
import csv
import gzip
from collections import defaultdict

target_species_ids = {"4911", "1003335", "1587", "405566", "585520", "767456", "267818", "190906", "33962", "4909"}
read_to_species = {}
for kraken_file in sys.argv[1:]:
    print(f"Reading: {kraken_file}")
    with open(kraken_file, "rt") as f:
        for line in f:
            parts = line.strip().split("\t", 3)
            species_id = parts[2]
            read_id = parts[1]
            if species_id in target_species_ids:
                read_to_species[read_id] = species_id

total_hits = defaultdict(list)
total_num = 0
for kraken_file in sys.argv[1:]:
    fastq_file = kraken_file.replace("barcode0", "barcode").replace(".kraken.core.out", ".fastq.gz")
    print(f"Reading: {fastq_file}")
    with gzip.open(fastq_file, "rt") as f:
        while True:
            header = f.readline().strip()
            seq = f.readline().strip()
            f.readline()  # plus line
            f.readline()  # quality line

            if not header:
                break

            if header.startswith("@"):
                read_id = header[1:].split()[0]
                species_id = read_to_species.get(read_id)
                if species_id is not None:
                    total_hits[species_id].append(f">{read_id}\n{seq}")

            total_num += 1
            if total_num % 200000 == 0:
                print(".")

for species_id, content in total_hits.items():
    filename = f"{species_id}.metagenomic.fasta"
    print(f"Species: {species_id}")
    print(f"Outputting to: {filename}")
    with open(filename, "w") as f:
        f.write("\n".join(content) + "\n")
