import subprocess
import re

def align(first, second):
    retaln = ""
    for a, b in zip(first, second):
        if a == '-' or b == '-':
            retaln += " "
        elif a == b:
            retaln += "|"
        else:
            retaln += " "
    return retaln

blast_cmd_template = (
    "blastn -query XXXX -db YYYY "
    "-outfmt '7 qseqid sseqid pident mismatch gapopen qstart qend "
    "sstart send evalue bitscore length qlen qseq sseq' "
    "-word_size 8 -evalue 1e-5 -num_threads 12"
)

with open("targets.txt", "r") as infile:
    for line in infile:
        line = line.strip()
        if not line:
            continue
        species, db = line.split("\t")
        species_file = species.replace(" ", "_") + ".metagenomic.fasta"
        print(f"{species_file} {db}")

        blast_cmd = blast_cmd_template.replace("XXXX", species_file).replace("YYYY", db)
        seen = set()
        total_hits = 0
        total_seqs = 0

        with open(f"{species_file}.blasthits.txt", "w") as outfile:
            process = subprocess.Popen(blast_cmd, shell=True, stdout=subprocess.PIPE, text=True)
            current_query = None
            for line in process.stdout:
                line = line.strip()
                if not line:
                    continue

                query_match = re.match(r"^# Query: (\S+)", line)
                if query_match:
                    current_query = query_match.group(1)
                    continue

                hits_match = re.match(r"^# (\d+) hits found", line)
                if hits_match:
                    if int(hits_match.group(1)) == 0:
                        outfile.write(f"{current_query}\tNo Hits\tNA\n")
                        total_seqs += 1
                    continue

                if not line.startswith("#") and current_query and current_query not in seen:
                    parts = line.split("\t")
                    qfrac = float(parts[11]) / float(parts[12])
                    qfrac_str = f"{qfrac * 100:.2f}%"
                    outfile.write(f"{current_query}\t{parts[1]}\t{parts[2]}\t{qfrac_str} [{parts[11]} / {parts[12]}]\t{parts[7]}\t{parts[8]}\t{parts[5]}\t{parts[6]}\n")

                    aln_marker = align(parts[13], parts[14])
                    qseq_chunks = [parts[13][i:i+80] for i in range(0, len(parts[13]), 80)]
                    sseq_chunks = [parts[14][i:i+80] for i in range(0, len(parts[14]), 80)]
                    aln_chunks = [aln_marker[i:i+80] for i in range(0, len(aln_marker), 80)]

                    for q, a, s in zip(qseq_chunks, aln_chunks, sseq_chunks):
                        outfile.write(f"#\n#[aln]\t{q}\n#[aln]\t{a}\n#[aln]\t{s}\n")
                    outfile.write("#\n")

                    seen.add(current_query)
                    total_hits += 1
                    total_seqs += 1

        print(f"Total Hits: {total_hits} in {total_seqs} (seqs)")
