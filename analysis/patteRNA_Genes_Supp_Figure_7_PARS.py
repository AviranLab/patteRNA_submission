import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# DATAPATH
pars_files = ["PARS/5SL_d75/scores.txt",
              "PARS/4SL_d75/scores.txt"]
pars_labels = ["5SL",
               "4SL"]

shape_files = ["genome_scores/scores_genome_1M6_full_length.txt",
               "genome_scores/scores_genome_1M7_full_length.txt",
               "genome_scores/scores_genome_NMIA_full_length.txt",
               "genome_scores/scores_genome_2009_full_length.txt"]
shape_labels = ["HIV 1M6",
                "HIV 1M7",
                "HIV NMIA",
                "HIV 2009"]
shape_id = 1
top_n = 1
linestyles = ['-', '--', ':']


def read_scores(fp):
    scores = []

    with open(fp, "r") as f:
        _ = f.readline()

        for line in f:

            fields = line.strip().split()
            curr_score = np.float(fields[4])

            if ~np.isnan(curr_score):
                scores.append(curr_score)

    scores = np.array(scores, dtype=float)

    return scores


def cdf(x):
    n = len(x)
    x = np.sort(x)
    y = np.array(range(n)) / float(n)

    return x, y


def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Read SHAPE
    shape = []
    shape_anchors = []

    for file in shape_files:
        curr_score = read_scores(file)
        shape.append(curr_score)
        shape_anchors.append(curr_score[0])

    # Read PARS
    pars = []
    for file in pars_files:
        curr_score = read_scores(file)
        pars.append(curr_score)

        for i, anchor in enumerate(shape_anchors):
            print(file, shape_labels[i], "{:d}/{:d}".format(np.sum(curr_score > anchor) + 1, len(curr_score)))

    # plot SHAPE
    x, y = cdf(shape[shape_id])
    y = 1 - y
    ax.plot(x, y,
            color="black",
            label=shape_labels[shape_id] + " (n = {:d})".format(len(shape[shape_id])))

    ax.add_patch(patches.Rectangle((10.6, -0.05),
                                   13.2-10.6,
                                   1.1,
                                   alpha=0.3,
                                   color="black",
                                   edgecolor=None))

    for i, score in enumerate(pars):
        x, y = cdf(pars[i])
        y = 1 - y
        plt.plot(x, y,
                 linestyle=linestyles[i],
                 color="grey",
                 label="{} RRE in PARS (n = {:d})".format(pars_labels[i], len(score)))

    plt.xlabel("c-score cutoff")
    plt.ylabel("1-CDF")
    plt.legend()
    plt.tight_layout()

    fig.savefig("figure.png", dpi=300)


if __name__ == "__main__":
    main()
