#!/usr/bin/env python3

from pathlib import Path
import sys

id = sys.argv[1]
genome_size_num = sys.argv[2]

stats_file = Path("nanoplot/NanoStats.txt")
stats = {}

for line in stats_file.read_text().splitlines():
    if ":" in line:
        key, value = line.split(":", 1)
        stats[key.strip()] = value.strip()
total_bases = float(stats["Number of bases"])
mean_q = float(stats["Mean read quality"])
read_n50 = float(stats["Read length N50"])

genome_size = genome_size_num
ont_depth = total_bases / genome_size

with open("${meta.id}_ont_quality.tsv", "w") as out:
    out.write("sample\\tont_depth\\tont_mean_q\\tont_n50\\n")
    out.write( f"${id}\\t{ont_depth:.2f}\\t{mean_q:.2f}\\t{read_n50:.0f}\\n")