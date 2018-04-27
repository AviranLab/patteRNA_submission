import sys
import os
import numpy as np
import copy as cp
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(".."))

from analysis import bundlelib

# patteRNA_process_motif_bundle_average_short.py
"""The following script produces useful plots when assessing the scores of multiple motifs (bundles)
against multiple transcripts. This version in particular focuses on short motif scores by averaging
over a bundle of scores for each start position."""

# score_file = "scores_no_bundles_no_log.txt"
score_file = "scores/scores_mixtures.txt"
ensemble_file = "motif_files/RRE_short.txt"

# Data structure used when transcripts are read into memory: dict with motifs as keys
# Handy map for motif <--> path
empty_transcript, motif_map, _ = bundlelib.read_ensemble(ensemble_file)

# Read input file
print("Compiling unsorted scores...")
data_inv = bundlelib.compile_bundle_dict(score_file, empty_transcript, motif_map, path_column=6)

print("... done.")

xmin = 155
x = 163
xmax = 170

diff = 60  # Hardcoded offset between data index and nucleotide position

# Each figure refers to a transcript
# Each subplot refers to a different motif
print("Preparing plots...")
mutants = {"4SL": None, "4SL90-5SL10": None, "4SL80-5SL20": None, "4SL70-5SL30": None, "4SL60-5SL40": None,
           "4SL50-5SL50": None, "4SL40-5SL60": None, "4SL30-5SL70": None,
            "4SL20-5SL80": None, "4SL10-5SL90": None, "5SL": None}

# mutants = {"4SL": None, "4SL80-5SL20": None, "4SL60-5SL40": None,
#            "4SL50-5SL50": None, "4SL40-5SL60": None,
#             "4SL20-5SL80": None, "5SL": None}

color = None
xdata = []
ydata = []
switched = None
for i, transcript in zip(range(len(mutants.keys())), sorted(mutants.keys())):

    mutants[transcript] = []

    for j, mutant in zip(range(2), {"mutA", "mutB"}):

        x = i

        if i == 0:
            x = 10

        if i == 10:
            x = 0

        if j == 0:
            shift = -0.4
            color = (0.3, 0.3, 0.3)
        if j == 1:
            shift = 0.4
            color = (0.7, 0.7, 0.7)

        mutants[transcript].append(np.mean(data_inv[transcript][mutant][163-60]))
        # plt.bar(2*x+shift, np.mean(data_inv[transcript][mutant][163-60]), fc=color, ec="black", label=mutant)
        # plt.errorbar(2*x+shift,
        #              np.mean(data_inv[transcript][mutant][163-60]),
        #              np.std(data_inv[transcript][mutant][163-60]),
        #              color="black")
    xdata.append(x)
    if switched is None:
        if mutants[transcript][0]/mutants[transcript][1] > 1:
            switched=False
        else:
            switched=True

    if switched:
        ydata.append(mutants[transcript][0]/mutants[transcript][1])
    else:
        ydata.append(mutants[transcript][1] / mutants[transcript][0])


plt.figure(1)
shift=0.04
xdata = np.multiply(2, xdata)
ind = np.argsort(xdata)
plt.plot([xdata[i] for i in ind], [ydata[i]+shift for i in ind], "s-k")
plt.xlabel("SHAPE Ratio")
plt.xticks((0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20), ("100% 5SL",
                                                     "90/10",
                                                     "80/20",
                                                     "70/30",
                                                     "60/40",
                                                     "50/50",
                                                     "40/60",
                                                     "30/70",
                                                     "20/80",
                                                     "10/90",
                                                     "100% 4SL"), rotation="45")
for tick in plt.gca().xaxis.get_majorticklabels():
    tick.set_horizontalalignment("right")

plt.ylabel("$\\frac{c_{\\mathrm{5SL}}}{c_{\\mathrm{4SL}}}$", rotation="horizontal", labelpad=30)
plt.gca().yaxis.label.set_size(20)
plt.axhline(1, color="black")
plt.axhline(0.47, color="black", ls='--')
plt.ylim((0, 2))
plt.yticks((0, 1, 2))
plt.text(0, 0.5, "$Genomic$ $Ratio$ (NMIA)")
plt.tight_layout()
plt.show()
