import sys
import csv
import gzip
from collections import defaultdict

# Define target species IDs
target_species_ids = {
    "267818", "190906"
}

# Mapping from read ID to species ID
read_to_species = {}

# First pass: read Kraken outputs
for kraken_file in sys.argv[1:]:
    print(f"Reading: {kraken_file}")
    with open(kraken_file, "rt") as f:
        for line in f:
            parts = line.strip().split("\t", 3)
            species_id = parts[2]
            read_id = parts[1]
            if species_id in target_species_ids:
                read_to_species[read_id] = species_id

# Second pass: match reads and collect FASTQ entries
total_hits = defaultdict(list)
total_num = 0

for kraken_file in sys.argv[1:]:
    fastq_file = kraken_file.replace("barcode0", "barcode").replace(".kraken.core.out", ".fastq.gz")
    print(f"Reading: {fastq_file}")
    with gzip.open(fastq_file, "rt") as f:
        while True:
            header = f.readline().strip()
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            if not header:
                break

            if header.startswith("@"):
                read_id = header[1:].split()[0]
                species_id = read_to_species.get(read_id)
                if species_id is not None:
                    fastq_entry = f"{header}\n{seq}\n{plus}\n{qual}"
                    total_hits[species_id].append(fastq_entry)

            total_num += 1
            if total_num % 200000 == 0:
                print(f"Processed {total_num} reads...")

# Write output FASTQ files
for species_id, content in total_hits.items():
    filename = f"{species_id}.metagenomic.fastq"
    print(f"Species: {species_id}")
    print(f"Outputting to: {filename}")
    with open(filename, "w") as f:
        f.write("\n".join(content) + "\n")