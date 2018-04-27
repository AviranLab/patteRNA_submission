import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta, norm, dgamma, exponnorm, pearson3, genextreme, genlogistic, kappa4, skewnorm
from copy import deepcopy
sys.path.insert(0, os.path.abspath(".."))
from analysis import bundlelib

null_file = "sample_nulls/h0.temp"
null_ensemble_file = "motif_files/normalization_motifs.txt"
score_file = "scores/scores_normalization.txt"

empty_transcript, motif_map, motif_map_r = bundlelib.read_ensemble(null_ensemble_file)

null_score_repo = bundlelib.read_nulls(null_file, deepcopy(empty_transcript), motif_map)

score_repo = bundlelib.compile_bundle_dict(score_file, deepcopy(empty_transcript), motif_map, score_column=3, transcript_wise=False)
c_score_repo = bundlelib.compile_bundle_dict(score_file, deepcopy(empty_transcript), motif_map, score_column=4, transcript_wise=False)


def plot_beta_fit(data, line_spec='-', x=np.linspace(-350, 350, 1000)):

    min_score = np.min(data)-1
    max_score = np.max(data)+1

    floc = min_score
    fscale = 2*(max_score-min_score)
    fscale+=fscale

    # a, b, loc, scale = beta.fit(data, floc=floc, fscale=fscale)
    # y = beta.pdf(x, a=a, b=b, loc=floc, scale=fscale)
    #
    # mu, std = norm.fit(data)
    # yn = norm.pdf(x, mu, std)
    #
    # K, loc, scale=pearson3.fit(data, loc=floc, scale=scale)
    # yd = pearson3.pdf(x, skew=K, loc=loc, scale=scale)
    #
    K, loc, scale = genlogistic.fit(data)
    yl = genlogistic.pdf(x, c=K, loc=loc, scale=scale)


    plt.plot(x, yl, line_spec, color='black', alpha=1)
    print(genlogistic.stats(c=K, moments='mvsk'))


bins = np.linspace(-200, 200, 50)
plt.figure(1)

plt.subplot(3, 3, 1)
plt.hist(null_score_repo["one"], bins, normed=1, color="grey")
plot_beta_fit(null_score_repo["one"])
plt.xlim((-150, 150))
plt.xticks((-100, 0, 100))
plt.ylim((0, 0.02))
plt.yticks((0, 0.02), ("0", "0.02"))
plt.title("$H_0$")
plt.ylabel("Density")
plt.gca().yaxis.set_label_coords(-0.1, 0.5)
plt.subplot(3, 3, 2)
plt.hist(score_repo['all']['one'], bins, normed=1, color="grey")
plot_beta_fit(null_score_repo["one"], '--')
plt.title("Putative Sites")
plt.xlim((-150, 150))
plt.xticks((-100, 0, 100))
plt.ylim((0, 0.02))
plt.yticks(())
plt.subplot(3, 3, 3)
plt.hist(c_score_repo['all']['one'], np.linspace(0, 4, 60), normed=1, color="red")
plt.ylim((0, 2))
plt.yticks((0, 2), ("0", "1"))
plt.ylabel("Density")
plt.xlim((0, 4))
plt.xticks((0, 1, 2, 3, 4))
plt.title("Normalized\nTarget Scores")

plt.subplot(3, 3, 4)
plt.hist(null_score_repo["two"], bins, normed=1, color="grey")
plot_beta_fit(null_score_repo["two"])
plt.xlim((-150, 150))
plt.xticks((-100, 0, 100))
plt.ylim((0, 0.02))
plt.ylabel("Density")
plt.yticks((0, 0.02), ("0", "0.02"))
plt.gca().yaxis.set_label_coords(-0.1, 0.5)
plt.subplot(3, 3, 5)
plt.hist(score_repo['all']['two'], bins, normed=1, color="grey")
plot_beta_fit(null_score_repo["two"], '--')
plt.xlim((-150, 150))
plt.xticks((-100, 0, 100))
plt.ylim((0, 0.02))
plt.yticks(())
plt.subplot(3, 3, 6)
plt.hist(c_score_repo['all']['two'], np.linspace(0, 4, 60), normed=1, color="red")
plt.ylim((0, 2))
plt.yticks((0, 2), ("0", "1"))
plt.ylabel("Density")
plt.xlim((0, 4))
plt.xticks((0, 1, 2, 3, 4))

plt.subplot(3, 3, 7)
plt.hist(null_score_repo["thr"], bins, normed=1, color="grey")
plot_beta_fit(null_score_repo["thr"])
plt.xlim((-150, 150))
plt.xticks((-100, 0, 100))
plt.ylim((0, 0.02))
plt.xlabel("Raw Score")
plt.ylabel("Density")
plt.yticks((0, 0.02), ("0", "0.02"))
plt.gca().yaxis.set_label_coords(-0.1, 0.5)
plt.subplot(3, 3, 8)
plt.hist(score_repo['all']['thr'], bins, normed=1, color="grey")
plot_beta_fit(null_score_repo["thr"], '--')
plt.xlim((-150, 150))
plt.xticks((-100, 0, 100))
plt.xlabel("Raw Score")
plt.ylim((0, 0.02))
plt.yticks(())
plt.subplot(3, 3, 9)
plt.hist(c_score_repo['all']['thr'], np.linspace(0, 4, 60), normed=1, color="red")
plt.ylim((0, 2))
plt.yticks((0, 2), ("0", "1"))
plt.ylabel("Density")
plt.xlim((0, 4))
plt.xticks((0, 1, 2, 3, 4))
plt.xlabel("c-score")


plt.subplots_adjust(right=0.8, wspace=0, hspace=0.25)
plt.subplot(3, 3, 3)
pos = plt.gca().get_position()
plt.gca().set_position([pos.x0+0.1, pos.y0, pos.width, pos.height])
plt.gca().yaxis.set_label_coords(-0.1, 0.5)
plt.subplot(3, 3, 6)
pos = plt.gca().get_position()
plt.gca().set_position([pos.x0+0.1, pos.y0, pos.width, pos.height])
plt.gca().yaxis.set_label_coords(-0.1, 0.5)
plt.subplot(3, 3, 9)
pos = plt.gca().get_position()
plt.gca().set_position([pos.x0+0.1, pos.y0, pos.width, pos.height])
plt.gca().yaxis.set_label_coords(-0.1, 0.5)

plt.show()
