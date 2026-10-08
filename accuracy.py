#!/usr/bin/env/python3

import sys

total_length=0
total_error_count = 0

with open(sys.argv[1], "rt") as f:
    for line in f:
        field=line.slipt("\t")
        percent_identity=float(field[2])
        alignment_length=int(field[3])
        if percent_identity < 95.0:
            continue
        elif alignment_length < 1000:
            continue
        else:
            total_length += alignment_length
            error_count=int(alignment_length * (100.0 - percent_identity) / 100.0)
            total_error_count += error_count
            complete_percent_identity=100.0 * (1.0-(total_error_count/total_length))
            complete_error_rate = 100.0 - complete_percent_identity
    print("%.2f" % complete_percent_identity + "\t" + "%.2f" % complete_error_rate)
