import sys
import os
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(".."))
from analysis import bundlelib

"""patteRNA_process_motif_bundle_average_full.py"""
"""The following script produces useful plots when assessing the scores of multiple motifs (bundles)
against multiple transcripts. This version uses the full-length RNA motif, and therefore only
one start position is possible."""

score_files = {"genome_scores/scores_genome_1M6_full_length.txt",
               "genome_scores/scores_genome_1M7_full_length.txt",
               "genome_scores/scores_genome_NMIA_full_length.txt",
               "genome_scores/scores_genome_2009_full_length.txt"}

# Data structure used when transcripts are read into memory: dict with motifs as keys
# Handy map for motif <--> path
empty_transcript, motif_map, motif_map_r = bundlelib.read_ensemble("motif_files/RRE_full.txt")
data = dict()
bars5 = dict()
bars4 = dict()

print("Reading score files...")
for score_file in score_files:

    orig_data = bundlelib.compile_bundle_dict(score_file, empty_transcript, motif_map, transcript_wise=1)

    data[score_file] = None
    bars5[score_file] = []
    bars4[score_file] = []
    for key in orig_data.keys():
        data[score_file] = orig_data[key]

    for score in data[score_file]["mutA"]:
        bars5[score_file].append(np.mean(score))

    for score in data[score_file]["mutB"]:
        bars4[score_file].append(np.mean(score))

print("... done.")
# Each figure refers to a transcript
# Each subplot refers to a different motif
print("Preparing plots...")


print(len(sorted(motif_map.keys())[0]))

plt.figure(num=None, figsize=(10, 10))
i = 1
j = 2
for score_file in score_files:

    print(score_file)
    plt.subplot(4, 2, i)
    i += 2
    plt.scatter(np.array(range(len(bars5[score_file])))+1, bars5[score_file],

                color="darkblue", s=np.power(np.add(bars5[score_file], 1), 1.1)-1)
    plt.xlim((1, 9174-len(sorted(motif_map.keys())[0])+2))
    plt.xticks((1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9174-len(sorted(motif_map.keys())[0])+2))
    plt.ylim((0, 15))
    plt.yticks((0, 5, 10, 15))
    plt.grid(axis='y')
    plt.ylabel("$c$-score")
    if i == 9:
        plt.xlabel("Position on HIV-1 Genome")

    plt.subplot(4, 2, j)
    j += 2
    plt.scatter(np.array(range(len(bars4[score_file]))) + 1, bars4[score_file],
                color="darkblue", s=np.power(np.add(bars4[score_file], 1), 1.1) - 1)
    plt.xlim((1, 9174 - len(sorted(motif_map.keys())[0])+2))
    plt.xticks((1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9174-len(sorted(motif_map.keys())[0])+2))
    plt.ylim((0, 15))
    plt.yticks((0, 5, 10, 15))
    plt.grid(axis='y')
    if j == 10:
        plt.xlabel("Position on HIV-1 Genome")

plt.tight_layout()
plt.show()