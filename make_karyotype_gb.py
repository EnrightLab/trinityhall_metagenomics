import random
import sys
import os
import re

def rand_col():
    return tuple(random.randint(0, 255) for _ in range(3))


def avg_array(coverage, a, b):
    total = sum(coverage[i] for i in range(a, b + 1))
    return total / (b - a + 1)


def parse_locus_size(filepath):
    with open(filepath) as f:
        for line in f:
            if line.startswith("LOCUS"):
                parts = line.split()
                return int(parts[2])
    return 0

def parse_location(loc_str):
    strand = "+"
    if loc_str.startswith("complement"):
        strand = "-"
        loc_str = re.sub(r"complement\((.*)\)", r"\1", loc_str)

    # Handle join(...) or order(...)
    if loc_str.startswith("join") or loc_str.startswith("order"):
        loc_str = re.sub(r"(join|order)\((.*)\)", r"\2", loc_str)

    segments = []
    for part in loc_str.split(","):
        if ".." in part:
            start, end = part.split("..")
            start = int(start.replace("<", "").replace(">", ""))
            end = int(end.replace("<", "").replace(">", ""))
            segments.append((start, end))
        else:
            pos = int(part.replace("<", "").replace(">", ""))
            segments.append((pos, pos))

    return segments, strand

def write_karyotype(size, chunksize, remainder, filekar):
    filekar.write(f"chr\t-\tchr1\t1\t0\t{size}\tblack\n")
    for i in range(size // chunksize):
        color = "black" if i % 2 else "grey"
        start = i * chunksize
        end = start + chunksize
        filekar.write(f"band\tchr1\t1.1\t1.1\t{start}\t{end}\t{color}\n")
    start = end
    end = start + remainder
    filekar.write(f"band\tchr1\t1.1\t1.1\t{start}\t{end}\tred\n")


def parse_hits(filepath, coverage):
    largest = 0
    max_val = 0
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            start, end = int(parts[8]), int(parts[9])
            if start > end:
                start, end = end, start
            largest = max(largest, end)
            for i in range(start, end):
                coverage[i] += 1
                max_val = max(max_val, coverage[i])
    return largest, max_val


def write_coverage_track(largest, small_chunk_size, coverage, max_val, size, filehits):
    avg_coverage = 0
    coverage_i = 0
    for i in range(largest // small_chunk_size):
        start = i * small_chunk_size
        end = start + small_chunk_size - 1
        avg = avg_array(coverage, start, end)
        filehits.write(f"chr1\t{start}\t{end}\t{avg / max_val}\n")
        avg_coverage += avg
        coverage_i += 1

    start = end + 1
    remainder = largest % small_chunk_size
    end = start + remainder - 1
    avg = avg_array(coverage, start, end)
    filehits.write(f"chr1\t{start}\t{end}\t{avg}\n")


def write_labels_and_highlights(genbank_file, filelab1, filelab2, fileg1, fileg2, all_names):
    with open(genbank_file) as f:
        in_feature_table = False
        label = ""
        segments = []
        strand = "+"
        for line in f:
            if line.startswith("FEATURES"):
                in_feature_table = True
                continue
            if not in_feature_table:
                continue
            if line.startswith("ORIGIN"):
                break
            if line.startswith("     gene"):
                tokens = line.split()
                if len(tokens) > 1:
                    loc = tokens[1]
                    try:
                        segments, strand = parse_location(loc)
                    except Exception as e:
                        print(f"Skipping unrecognized gene location: {loc} ({e})")
                        continue
            elif "/gene=" in line or "/product=" in line:
                label = line.split("=")[1].strip().strip('"')
                if label and segments:
                    labelfile = filelab1 if strand == "+" else filelab2
                    highfile = fileg1 if strand == "+" else fileg2
                    for start, end in segments:
                        labelfile.write(f"chr1\t{start}\t{end}\t{label}\n")
                        if all_names:
                            highfile.write(f"chr1\t{start}\t{end}\tfill_color=red\n")
                    segments = []

def main(genbank_file, hits_file, skew_file=None, all_names=False):
    size = parse_locus_size(genbank_file)
    ideal = size // 30
    chunksize = ideal if size < 1_000_000 else 100_000
    small_chunk_size = max(ideal // 100, 1000)
    chunks = size // chunksize
    remainder = size % chunksize

    with open("karyotype.txt", "w") as filekar, \
         open("data_tile.txt", "w") as filehits, \
         open("labels_forward.txt", "w") as filelab1, \
         open("labels_reverse.txt", "w") as filelab2, \
         open("highlights_forward.txt", "w") as fileg1, \
         open("highlights_reverse.txt", "w") as fileg2:

        write_karyotype(size, chunksize, remainder, filekar)
        write_labels_and_highlights(genbank_file, filelab1, filelab2, fileg1, fileg2, all_names)

        coverage = [0] * (size + 1)
        largest, max_val = parse_hits(hits_file, coverage)
        write_coverage_track(largest, small_chunk_size, coverage, max_val, size, filehits)

    if skew_file and os.path.exists(skew_file):
        with open(skew_file) as f, open("data_skew.txt", "w") as fileskew:
            for line in f:
                if line.startswith("Sequence"):
                    continue
                parts = line.strip().split("\t")
                pos = int(parts[1])
                val = parts[2]
                fileskew.write(f"chr1\t{pos}\t{pos + 20000 - 1}\t{val}\n")
            fileskew.write(f"chr1\t{pos + 20000}\t{size}\t{val}\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <genbank_file> <hits_file> [skew_file] [--names]")
        sys.exit(1)

    genbank_file = sys.argv[1]
    hits_file = sys.argv[2]

    skew_file = None
    all_names = False
    args = sys.argv[3:]

    for arg in args:
        if arg == "--names":
            all_names = True
        elif not arg.startswith("--"):
            skew_file = arg

    main(genbank_file, hits_file, skew_file, all_names)