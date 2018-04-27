import sys
import os
import numpy as np
import copy as cp
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(".."))

from analysis import bundlelib

"""patteRNA_process_motif_bundle_average_full.py"""
"""The following script produces useful plots when assessing the scores of multiple motifs (bundles)
against multiple transcripts. This version uses the full-length RNA motif, and therefore only
one start position is possible."""

score_file = "genome_scores/scores_genome_NMIA_full_length.txt"

# Data structure used when transcripts are read into memory: dict with motifs as keys
# Handy map for motif <--> path
empty_transcript, motif_map, motif_map_r = bundlelib.read_ensemble("motif_files/RRE_full.txt")
orig_data = bundlelib.compile_bundle_dict(score_file, empty_transcript, motif_map, transcript_wise=1)
# print("... done.")

data = {"HIV": None}
diff = 0  # Hardcoded offset between data index and nucleotide position
for key in orig_data.keys():
    data["HIV"] = orig_data[key]


# Each figure refers to a transcript
# Each subplot refers to a different motif
print("Preparing plots...")

bars = []
for score in data["HIV"]["mutB"]:
    bars.append(np.mean(score))

bars2 = []
for score in data["HIV"]["mutA"]:
    bars2.append(np.mean(score))

plt.figure(num=None, figsize=(7, 3))
plt.subplot(1, 3, (1, 2))
plt.scatter(np.array(range(len(bars)))+1, bars, color="darkblue", s=np.power(np.add(bars, 1), 1.1)-1)
#plt.scatter(np.array(range(len(bars)))+1, bars, color="darkblue", s=2)
plt.xlim((1, 9174-len(sorted(motif_map.keys())[0])))
plt.ylabel("$c$-score")
plt.grid(axis='y')
plt.xticks((1, 7306, 9174-len(sorted(motif_map.keys())[0])), ("1", "7306", str(2+9174-len(sorted(motif_map.keys())[0]))))
plt.xlabel("Position on HIV Genome")
plt.ylim((0, 15))
plt.yticks((0, 5, 10, 15))
#plt.text(184, bars2[183], "DIS")
plt.annotate("DIS", xy=(170, bars[168]+0.5), xytext=(183+50, bars[168]+2),
             arrowprops=dict(arrowstyle='->'))

#plt.annotate("$gag$-$pol$ Frameshift", xy=(1632, bars[1631]+0.5), xytext=(500, bars[1630]+1.6),
#             arrowprops=dict(arrowstyle='->'), size=7)

plt.annotate("RT$_{PK}$", xy=(2633, bars[2632]+0.5), xytext=(2633+50, bars[2632]+2),
             arrowprops=dict(arrowstyle='->'))

plt.annotate("ESSV\nJunction", xy=(4792, bars[4791]+0.1), xytext=(4792-250, bars[4791]+2),
             arrowprops=dict(arrowstyle='->'), size=9, horizontalalignment='center')

plt.annotate("RRE", xy=(7305-50, bars[7305]+0.1), xytext=(7305-1500, bars[7305]),
             arrowprops=dict(arrowstyle='->'))

plt.annotate("3$'$-TAR", xy=(8886, bars[8885]+0.1), xytext=(8886-1100, bars[8885]+2),
             arrowprops=dict(arrowstyle='->'))

plt.subplot(1, 3, 3)
plt.plot(np.array(range(len(bars)))+1, bars, color="darkblue", lw=1.5)
plt.xlim((7200, 7400))
plt.ylabel("$c$-score")
plt.grid(axis='y')
plt.xticks((7200, 7306, 7400))
plt.xlabel("Position on HIV Genome")
plt.ylim((0, 15))
plt.yticks((0, 5, 10, 15))

plt.tight_layout()
print("... done.")
plt.show()

