import sys
import os
import numpy as np
import copy as cp
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(".."))

from analysis import bundlelib
plt.rc('legend', fontsize=9)  # using a size in points


# patteRNA_process_motif_bundle_average_short.py
"""The following script produces useful plots when assessing the scores of multiple motifs (bundles)
against multiple transcripts. This version in particular focuses on short motif scores by averaging
over a bundle of scores for each start position."""

score_file = "scores/scores_sherpa_set_full_length.txt"
ensemble_file = "motif_files/RRE_full.txt"

# Data structure used when transcripts are read into memory: dict with motifs as keys
# Handy map for motif <--> path
empty_transcript, motif_map, _ = bundlelib.read_ensemble(ensemble_file)

# Read input file
print("Compiling unsorted scores...")
data_inv = bundlelib.compile_bundle_dict(score_file, empty_transcript, motif_map, score_column=4)
print(data_inv)
print("... done.")

x = 0

diff = 0  # Hardcoded offset between data index and nucleotide position

# Each figure refers to a transcript
# Each subplot refers to a different motif
print("Preparing plots...")
mutants = {"mutA": None, "mutB": None, "mutC": None, "mutD": None, "mutE": None}
plt.figure(1, figsize=(6, 6))
for i, transcript in zip(range(5), sorted(mutants.keys())):

    if i == 0:
        plt.subplot2grid((3, 3), (0, 0))

    if i == 1:
        plt.subplot2grid((3, 3), (0, 1))

    if i == 2:
        plt.subplot2grid((3, 3), (0, 2))

    if i == 3:
        plt.subplot2grid((3, 3), (1, 0))

    if i == 4:
        plt.subplot2grid((3, 3), (2, 0))


    x_data = []

    for j, mutant in zip(range(7), sorted(empty_transcript.keys())):
        try:
            temp = np.mean(data_inv[transcript][mutant][x - diff])
        except IndexError:
            temp = 0

        x_data.append(temp)

    plt.bar(range(5), x_data, color="grey", ec="black")
    plt.text(i, x_data[i]+0.075, "$\\bigstar$", horizontalalignment="center")
    plt.xlim((-1, 5))

    if i == 4:
        plt.xlabel("Target Path")

    plt.ylim((1, 4))
    if i == 3:
        plt.ylabel("$c$-score")

    plt.xticks((0, 1, 2, 3, 4), ("A", "B", "C", "D", "E"), rotation="horizontal")
    #plt.axhline(0, color='black', lw=1)
    plt.title("Mutant "+transcript[3])
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

ax = plt.subplot2grid((3, 3), (1, 1), colspan=2, rowspan=2)
x_map = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
x_4SL = []
x_5SL = []
for i, mutant in zip(range(5), sorted(empty_transcript.keys())):
    try:
        temp1 = np.mean(data_inv["4SL"][mutant][x - diff])
        temp2 = np.mean(data_inv["5SL"][mutant][x - diff])
    except IndexError:
        temp1 = 0
        temp2 = 0
    x_4SL.append(temp1)
    x_5SL.append(temp2)

plt.bar(x_map-0.2, x_5SL, width=0.4, ec="black", fc="black")
plt.bar(x_map+0.2, x_4SL, width=0.4, ec="black", fc=(0.9, 0.9, 0.9))
plt.text(x_map[0]-0.2, x_5SL[0]+0.075, "$\\bigstar$", horizontalalignment="center")
plt.text(x_map[1]+0.2, x_4SL[1]+0.075, "$\\bigstar$", horizontalalignment="center")
plt.xticks((0, 1, 2, 3, 4), (" A (5SL)", " B (4SL)", "C", "D", "E"), rotation="horizontal")
plt.axhline(0, color='black', lw=1, label="_nolegend_")
plt.xlim((-1, 5))
#plt.ylabel("$c$-score")
plt.ylim((1, 4))
plt.yticks((1, 2, 3, 4))
plt.xlabel("Target Path")
plt.title("Native Isomers")
plt.legend(("5SL", "4SL"), fancybox=False)
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)

plt.tight_layout()
plt.subplots_adjust(hspace=0.6, wspace=0.25)
plt.show()

print("... done.")
