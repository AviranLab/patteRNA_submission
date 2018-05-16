import numpy as np
import matplotlib.pyplot as plt
import pickle
import warnings
from copy import deepcopy

# MATPLOTLIB CONFIGS
plt.rcParams.update({'font.size': 13,
                     'axes.linewidth': 1.2,
                     'xtick.major.width': 1.2,
                     'ytick.major.width': 1.2})
DPI = 300
LW = 4
COLORS = ["red", "blue"]  # Colors for [paired, unpaired] states

# INPUTS
in_file = "KL_div/hiv_1M7_partitioned.shape"
full_model = "KL_div/full_set/trained_model.pickle"
partial_model = "KL_div/partial_set/trained_model.pickle"
ref_id_fp = "KL_div/partial_set/ref_id.txt"


# CLASSES AND FUNCTIONS
class GMMHMM:
    def __init__(self):
        pass

    def load(self, fp):
        """Load a pickled trained model."""

        with open(fp, "rb") as f:
            pickled_dict = pickle.load(f)

        # Assign loaded attributes
        kwargs2attr_deep(self, pickled_dict)


def read_fastalike(fp):
    """Reads a fasta-like file and returns a dictionary."""

    rnas = {}
    tr_name = None
    content = ""

    with open(fp, "r") as f:

        while True:

            line = f.readline().strip()

            if not line:
                break

            if line.startswith(">"):
                if tr_name is not None:  # Store this transcript
                    rnas[tr_name] = content

                tr_name = line.split(">")[1].strip()  # Get the new transcript
                content = ""
            else:
                content += line

        # Append the last entry of the file
        if tr_name is not None:
            rnas[tr_name] = content

    return rnas


def read_observations(fp):
    """Reads a fasta-like formatted file of observations and returns a dictionary."""

    rnas = {}

    file_content = read_fastalike(fp)

    for tr_name, field in file_content.items():
        obs = field.strip().split()
        obs = [word.replace('NA', 'nan') for word in obs]  # Handle NA
        obs = np.array(obs, dtype=float)
        obs[np.isinf(obs)] = np.nan  # Handle infinite values
        rnas[tr_name] = obs

    return rnas


def elog(x):
    """"Smart log implementation.

    Returns -inf for all values <= 0.

    Args:
        x (np.array): input vector.

    Returns:
        y (np.array): log transformed vector.

    """

    x = np.array(x, dtype=float)

    if x.size == 1:
        if x <= 0:
            y = -np.inf
        else:
            y = np.log(x)
    else:
        valid_mask = x > 0
        y = np.empty_like(x)
        y[:] = -np.inf
        y[valid_mask] = np.log(x[valid_mask])

    return y


def kwargs2attr_deep(obj, kwargs):
    """Deep copy attributes values based on keyword arguments."""
    if kwargs is None:
        return

    try:
        assert isinstance(kwargs, dict)
    except AssertionError:
        print("Attempted to set arguments not using a dictionary.")
        return

    for key in kwargs.keys():
        setattr(obj, key, deepcopy(kwargs[key]))


def wnormpdf(x, mean=0, var=1, w=1, w_min=0):
    """Weighted Normal PDF.

    Add a small value (1e-20) if likelihood=0 because a PDF should never generate a "true" 0.

    Args:
        x: Vector of input values to compute the density at.
        mean: Mean.
        var: Variance.
        w: Weight.
        w_min: Minimum weight allowed.

    Returns:
        y: Likelihoods.
        w: Updated weight.

    """

    stdev = np.sqrt(var)

    # Set w to 0 if either sigma or w reach threshold values.
    if (stdev < 0) | (w <= w_min):
        y = x
        w = 0
    else:
        # noinspection PyTypeChecker
        u = np.array((x - mean) / stdev, dtype=float)
        y = np.exp(-u * u / 2) / (np.sqrt(2 * np.pi) * stdev)
        # noinspection PyTypeChecker
        if np.any(y == 0):
            y += 1e-20

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y *= w

    return y, w


def concatenate_across_rna(in_data, ref_id=None):
    if ref_id is None:
        ref_id = in_data.keys()

    v = []
    for k in ref_id:
        v += list(in_data[k])
    v = np.array(v)
    return v


def plot_fit(model, data, ax, ax_id, xlim=None, ylim=None):
    n_bins = 100

    obs = data
    obs = obs[(~np.isnan(obs)) & (~np.isinf(obs)) & (obs != 0)]

    # Get the probability of each states
    p_st = model.gmm_gamma
    y_all = None

    # Plot the entire data histogram
    if xlim is None:
        hist, bins = np.histogram(obs, bins=n_bins, density=True)
    else:
        hist, bins = np.histogram(obs, bins=n_bins, density=True, range=xlim)

    if ylim is not None:
        ax[ax_id].set_ylim(ylim)

    # Add the Gaussian component subplots for each state

    for i in model.states:

        y = None
        for m in range(model.K):
            y_sub, _ = wnormpdf(bins,
                                model.mu[i, m],
                                model.sigma[i, m],
                                model.w[i, m], 0)
            y_sub = y_sub * p_st[i]
            ax[ax_id].plot(bins, y_sub, color=COLORS[i], linestyle="--")
            y = y_sub if y is None else y + y_sub
        if i == 0:
            ax[ax_id].plot(bins, y, color=COLORS[i], linewidth=LW, label="unpaired")
        else:
            ax[ax_id].plot(bins, y, color=COLORS[i], linewidth=LW, label="paired")
        y_all = y if y_all is None else y_all + y

    # Whole data fit
    _, _, _ = ax[ax_id].hist(obs,
                             normed=True,
                             bins=bins,
                             label="data",
                             color="Grey")
    ax[ax_id].plot(bins, y_all, "k", linewidth=LW, label="patteRNA")


# PLOT
# Generate the figure
# noinspection PyTypeChecker
fig, axes = plt.subplots(nrows=1, ncols=2,
                         sharex=False,
                         sharey=False,
                         figsize=(12, 6),
                         dpi=DPI)
ax = axes

# FULL SET
shape = read_observations(in_file)
shape = concatenate_across_rna(shape)
shape = np.array(shape)
shape = elog(shape)

model = GMMHMM()
model.load(full_model)
# noinspection PyTypeChecker
plot_fit(model, shape, ax, ax_id=0, xlim=[-8, 2], ylim=[0, 0.5])

# PARTIAL SET
ref_id = []
with open(ref_id_fp, "r") as f:
    for line in f:
        line = line.strip()
        ref_id.append(line)

shape = read_observations(in_file)
shape = concatenate_across_rna(shape, ref_id=ref_id)
shape = np.array(shape)
shape = elog(shape)

model = GMMHMM()
model.load(partial_model)
# noinspection PyTypeChecker
plot_fit(model, shape, ax, ax_id=1, xlim=[-8, 2], ylim=[0, 0.5])

# AXIS AND LEGEND
ax[0].set_xlabel("SHAPE (log-transformed)")
ax[1].set_xlabel("SHAPE (log-transformed)")
ax[0].set_ylabel("Probability Density")

ax[0].set_title("Entire set (n = 9150)")
ax[1].set_title("Partial set (n = 4949, KL-div = $8.67 \cdot 10^{-3}$)")

handles, labels = ax[0].get_legend_handles_labels()
ax[0].legend(handles, labels, loc="upper left")

plt.tight_layout()
plt.savefig("figure.png")
plt.clf()
