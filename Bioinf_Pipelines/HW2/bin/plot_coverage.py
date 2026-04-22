#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--sample", required=True)
args = parser.parse_args()

positions = []
depths = []

with open(args.input) as f:
    for line in f:
        chrom, pos, depth = line.strip().split("\t")
        positions.append(int(pos))
        depths.append(int(depth))

plt.figure(figsize=(12, 4))
plt.plot(positions, depths, linewidth=0.8)
plt.xlabel("Position")
plt.ylabel("Coverage")
plt.title(f"Coverage plot: {args.sample}")
plt.tight_layout()
plt.savefig(args.output, dpi=200)