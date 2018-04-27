import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
sys.path.insert(0, os.path.abspath("../"))
from analysis import bundlelib  # noqa

def create_benchmark_figure(score_file, ensemble_file):

    empty_transcript, motif_map, motif_map_r = bundlelib.read_ensemble(ensemble_file)
    score_repo = bundlelib.compile_bundle_dict(score_file,
                                               deepcopy(empty_transcript),
                                               motif_map,
                                               transcript_wise=False,
                                               score_column=3)

    bins = np.linspace(-600, 500, 80)
    plt.rc('font', size=8)

    plt.figure()
    gs = plt.GridSpec(3, 4)
    plt.subplot(gs[0, 0:2])
    plt.title("Kernel: 01111111000011111110")
    plt.hist(score_repo['all']['1d'], bins, color=(0.1, 0.4, 0.0, 0.5))
    plt.hist(score_repo['all']['2d'], bins, color=(0.1, 0.6, 0.0, 0.5))
    plt.hist(score_repo['all']['3d'], bins, color=(0.1, 0.8, 0.0, 0.5))
    plt.hist(score_repo['all']['4d'], bins, color=(0.1, 1.0, 0.0, 0.5))
    plt.legend(sorted({"1x", "2x", "3x", "4x"}), fancybox=False)
    plt.ylabel("Counts")
    plt.ylim((0, 2000))
    plt.yticks((0, 2000), ("0", "2000"))
    plt.xlim((-500, 500))
    plt.gca().yaxis.set_label_coords(-0.05, 0.5)

    plt.subplot(gs[1, 0:2])
    plt.title("Kernel: 00011111000011111000")
    plt.hist(score_repo['all']['1f'], bins, color=(0.05, 0.05, 0.4, 0.5))
    plt.hist(score_repo['all']['2f'], bins, color=(0.05, 0.05, 0.6, 0.5))
    plt.hist(score_repo['all']['3f'], bins, color=(0.05, 0.05, 0.8, 0.5))
    plt.hist(score_repo['all']['4f'], bins, color=(0.05, 0.05, 1.0, 0.5))
    plt.legend(sorted({"1x", "2x", "3x", "4x"}), fancybox=False)
    plt.ylabel("Counts")
    plt.ylim((0, 2000))
    plt.yticks((0, 2000), ("0", "2000"))
    plt.xlim((-500, 500))
    plt.gca().yaxis.set_label_coords(-0.05, 0.5)

    plt.subplot(gs[2, 0:2])
    plt.title("Kernel: 00000111000011100000")
    plt.hist(score_repo['all']['1h'], bins, color=(0.4, 0.1, 0.0, 0.5))
    plt.hist(score_repo['all']['2h'], bins, color=(0.6, 0.1, 0.0, 0.5))
    plt.hist(score_repo['all']['3h'], bins, color=(0.8, 0.1, 0.0, 0.5))
    plt.hist(score_repo['all']['4h'], bins, color=(1.0, 0.1, 0.0, 0.5))
    plt.legend(sorted({"1x", "2x", "3x", "4x"}), fancybox=False)
    plt.xlabel("Raw Score")
    plt.ylabel("Counts")
    plt.ylim((0, 2000))
    plt.yticks((0, 2000), ("0", "2000"))
    plt.xlim((-500, 500))
    plt.gca().yaxis.set_label_coords(-0.05, 0.5)

    ax = plt.subplot(gs[0:, 2:])
    x = [20, 40, 60, 80]
    short_hairpins = [np.mean(score_repo['all']['1f']),
                      np.mean(score_repo['all']['2f']),
                      np.mean(score_repo['all']['3f']),
                      np.mean(score_repo['all']['4f'])]

    long_hairpins = [np.mean(score_repo['all']['1a']),
                     np.mean(score_repo['all']['1b']),
                     np.mean(score_repo['all']['1c']),
                     np.mean(score_repo['all']['1d'])]

    all_unpaired = [np.mean(score_repo['all']['1k']),
                    np.mean(score_repo['all']['2k']),
                    np.mean(score_repo['all']['3k']),
                    np.mean(score_repo['all']['4k'])]

    mot1 = [np.mean(score_repo['all']['1a']),
            np.mean(score_repo['all']['2a']),
            np.mean(score_repo['all']['3a']),
            np.mean(score_repo['all']['4a'])]

    mot2 = [np.mean(score_repo['all']['1b']),
            np.mean(score_repo['all']['2b']),
            np.mean(score_repo['all']['3b']),
            np.mean(score_repo['all']['4b'])]

    mot3 = [np.mean(score_repo['all']['1c']),
            np.mean(score_repo['all']['2c']),
            np.mean(score_repo['all']['3c']),
            np.mean(score_repo['all']['4c'])]

    mot4 = [np.mean(score_repo['all']['1d']),
            np.mean(score_repo['all']['2d']),
            np.mean(score_repo['all']['3d']),
            np.mean(score_repo['all']['4d'])]

    mot5 = [np.mean(score_repo['all']['1e']),
            np.mean(score_repo['all']['2e']),
            np.mean(score_repo['all']['3e']),
            np.mean(score_repo['all']['4e'])]

    mot6 = [np.mean(score_repo['all']['1f']),
            np.mean(score_repo['all']['2f']),
            np.mean(score_repo['all']['3f']),
            np.mean(score_repo['all']['4f'])]

    mot7 = [np.mean(score_repo['all']['1g']),
            np.mean(score_repo['all']['2g']),
            np.mean(score_repo['all']['3g']),
            np.mean(score_repo['all']['4g'])]

    mot8 = [np.mean(score_repo['all']['1h']),
            np.mean(score_repo['all']['2h']),
            np.mean(score_repo['all']['3h']),
            np.mean(score_repo['all']['4h'])]

    mot9 = [np.mean(score_repo['all']['1i']),
            np.mean(score_repo['all']['2i']),
            np.mean(score_repo['all']['3i']),
            np.mean(score_repo['all']['4i'])]

    mot10 = [np.mean(score_repo['all']['1j']),
             np.mean(score_repo['all']['2j']),
             np.mean(score_repo['all']['3j']),
             np.mean(score_repo['all']['4j'])]

    mot11 = [np.mean(score_repo['all']['1k']),
             np.mean(score_repo['all']['2k']),
             np.mean(score_repo['all']['3k']),
             np.mean(score_repo['all']['4k'])]

    plt.plot(x, mot1, '.', color='black')
    plt.plot(x, mot2, '.', color='black')
    plt.plot(x, mot3, '.', color='black')
    plt.plot(x, mot4, '.', color='green')
    plt.plot(x, mot5, '.', color='black')
    plt.plot(x, mot6, '.', color='blue')
    plt.plot(x, mot7, '.', color='black')
    plt.plot(x, mot8, '.', color='red')
    plt.plot(x, mot9, '.', color='black')
    plt.plot(x, mot10, '.', color='black')
    plt.plot(x, mot11, '.', color='black')

    m, b = np.polyfit(x, mot1, 1)
    plt.plot(x, np.multiply(m, x)+b, '--', color='black')
    plt.text(83, mot1[3]-4, "P:UP = 1")
    m, b = np.polyfit(x, mot2, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='black')
    plt.text(83, mot2[3]-4, "P:UP = 0.9")
    m, b = np.polyfit(x, mot3, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='black')
    plt.text(83, mot3[3]-4, "P:UP = 0.8")
    m, b = np.polyfit(x, mot4, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='green')
    plt.text(83, mot4[3]-4, "P:UP = 0.7")
    m, b = np.polyfit(x, mot5, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='black')
    plt.text(83, mot5[3]-4, "P:UP = 0.6")
    m, b = np.polyfit(x, mot6, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='blue')
    plt.text(83, mot6[3]-4, "P:UP = 0.5")
    m, b = np.polyfit(x, mot7, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='black')
    plt.text(83, mot7[3]-4, "P:UP = 0.4")
    m, b = np.polyfit(x, mot8, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='red')
    plt.text(83, mot8[3]-4, "P:UP = 0.3")
    m, b = np.polyfit(x, mot9, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='black')
    plt.text(83, mot9[3]-4, "P:UP = 0.2")
    m, b = np.polyfit(x, mot10, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='black')
    plt.text(83, mot10[3]-4, "P:UP = 0.1")
    m, b = np.polyfit(x, mot11, 1)
    plt.plot(x, np.multiply(m, x) + b, '--', color='black')
    plt.text(83, mot11[3]-4, "P:UP = 0")


    plt.xlabel("Motif Length")
    plt.ylabel("Mean Raw Score")
    plt.title("Target Composition Influence on Mean Scores")

    plt.ylim((-250, 250))
    plt.yticks((-200, -100, 0, 100, 200))
    plt.xlim((15, 110))
    plt.xticks((20, 40, 60, 80))
    plt.gca().yaxis.set_label_coords(-0.09, 0.5)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    args = sys.argv
    score_file = "scores/scores_ray_benchmark.txt"
    ensemble_file = "motif_files/benchmarking.txt"

    create_benchmark_figure(score_file, ensemble_file)
