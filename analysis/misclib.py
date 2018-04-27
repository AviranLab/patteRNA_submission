import math
from scipy import stats
import os
import subprocess
import time
from datetime import timedelta
import numpy as np
from copy import deepcopy


def ct2pair(d, context):
    n = max(d.keys())

    pair = {}
    for key, val in d.items():
        if context == 2:  # paired/unpaired split
            if val == -1:
                pair[key] = 0  # unpaired
            else:
                pair[key] = 1  # paired
        elif context == 3:  # helix-end/stacked/unpaired split
            if val == -1:
                pair[key] = 0  # unpaired
            elif (key == 0) or (key == n):
                pair[key] = 1  # helix-end paired
            elif (d[key - 1] != -1) and (
                        (d[key + 1]) != -1) and (
                        d[key - 1] == (d[key] + 1)) and (
                        d[key + 1] == (d[key] - 1)):
                pair[key] = 2  # stacked paired
            else:
                pair[key] = 1  # helix-end paired

    return dict2numpy(pair)  # this ensures to return a properly sorted array


def dict2numpy(d):
    n = len(d)
    arr = np.zeros(n)
    arr[:] = np.nan  # initialize the array to NaNs
    for k, v in d.items():
        arr[k] = v

    return arr


def read_fasta(fp):
    seq = None
    with open(fp, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line.startswith(">"):  # skip comment line
                line = list(line)
                seq = line

    return seq


def read_shape(fp):
    shape = np.loadtxt(fp)
    shape = shape[:, 1]
    shape[shape == -999] = np.nan  # replace -999 with NaNs

    return shape


def read_ct(fp, context=2):
    EOS = False
    with open(fp) as f:
        n = int(f.readline().strip().split()[0])  # get the RNA length from the header
        rows = (line.split() for line in f)
        seq = []  # initialize sequence
        pair = {}  # initialize pair information
        for row in rows:
            if int(row[0]) == n:
                if EOS:
                    break
                else:
                    EOS = True  # End of structure signal, must break next iteration
            seq.append(row[1])
            pair[int(row[2])] = int(row[4]) - 1  # -1 ensures motif_position starts at 0 (unpaired bases get -1)

    # Convert .ct pairing information into proper pairing given the structural context
    pair_out = ct2pair(pair, context)

    # Convert .ct pairing information into dot-bracket RNA structure notation
    dot_out = ct2dot(pair)

    return pair_out, dot_out, seq


def ct2dot(d):
    d_dot = {}
    for key, val in d.items():
        if val == -1:
            if key not in d_dot:
                d_dot[key] = "."
        elif key not in d_dot:
            d_dot[key] = "("
            d_dot[val] = ")"

    # Convert the dict to a str + ensure order is maintained
    dot = ""
    for key in sorted(d_dot):
        dot = dot + d_dot[key]

    return dot


def fold(in_fasta, ref_ct, out_ct, out_score, env):
    subprocess.call("echo $DATAPATH", env=env, shell=True)

    # Set commands
    fold_cmd = "RNAprob {} {}" \
        .format(in_fasta,
                out_ct)

    score_cmd = "scorer {} {} {}" \
        .format(out_ct,
                ref_ct,
                out_score)

    # Run commands
    print(fold_cmd)
    print(score_cmd)
    p1 = subprocess.Popen(fold_cmd, env=env, shell=True)
    p1.wait()

    p2 = subprocess.Popen(score_cmd, env=env, shell=True)
    p2.wait()


def fold_only(in_fasta, out_ct, env):
    subprocess.call("echo $DATAPATH", env=env, shell=True)

    # Set commands
    fold_cmd = "RNAprob {} {}" \
        .format(in_fasta,
                out_ct)

    # Run commands
    print(fold_cmd)
    p1 = subprocess.Popen(fold_cmd, env=env, shell=True)
    p1.wait()


def fold_w_constraints(in_fasta, in_shape, ref_ct, out_ct, out_score, env):
    subprocess.call("echo $DATAPATH", env=env, shell=True)

    # Set commands
    fold_cmd = "RNAprob {} {} -sh {} -2s" \
        .format(in_fasta,
                out_ct,
                in_shape)

    score_cmd = "scorer {} {} {}" \
        .format(out_ct,
                ref_ct,
                out_score)

    # Run commands
    print(fold_cmd)
    print(score_cmd)
    p1 = subprocess.Popen(fold_cmd, env=env, shell=True)
    p1.wait()

    p2 = subprocess.Popen(score_cmd, env=env, shell=True)
    p2.wait()


def fold_only_w_constraints(in_fasta, in_shape, out_ct, env):
    subprocess.call("echo $DATAPATH", env=env, shell=True)

    # Set commands
    fold_cmd = "RNAprob {} {} -sh {} -2s" \
        .format(in_fasta,
                out_ct,
                in_shape)

    # Run commands
    print(fold_cmd)
    p1 = subprocess.Popen(fold_cmd, env=env, shell=True)
    p1.wait()


def read_fold_score(in_):
    """Read a .score file from an accuracy assessment after RNA secondary structure prediction."""

    with open(in_, "r") as f:
        # skip headers
        for i in range(3):
            next(f)

        line = f.readline()
        line = line.split("=")
        sensitivity = float(line[1].split("%")[0].strip())

        line = f.readline()
        line = line.split("=")
        ppv = float(line[1].split("%")[0].strip())

        line = f.readline()
        line = line.split("=")
        mcc = float(line[1].split("%")[0].strip())

        next(f)
        next(f)

        line = f.readline()
        line = line.split("=")
        tp = int(line[1].strip())

        line = f.readline()
        line = line.split("=")
        fp = int(line[1].strip())

        line = f.readline()
        line = line.split("=")
        tn = int(line[1].strip())

        line = f.readline()
        line = line.split("=")
        fn = int(line[1].strip())

    return sensitivity, ppv, mcc, tp, fp, tn, fn


def combine_scores(score_dir, suffix, rna_names, out_file):
    sum_sensitivity = 0
    sum_ppv = 0
    sum_mcc = 0
    sum_tp = 0
    sum_fp = 0
    sum_tn = 0
    sum_fn = 0

    with open(out_file, "w") as f:
        for rna in rna_names:
            sensitivity, ppv, mcc, tp, fp, tn, fn = read_fold_score(os.path.join(score_dir, rna + suffix))

            f.write("{} {:.1f}\n".format(rna, mcc))

            sum_sensitivity += sensitivity
            sum_ppv += ppv
            sum_mcc += mcc
            sum_tp += tp
            sum_fp += fp
            sum_tn += tn
            sum_fn += fn

        # avg_sensitivity = sum_sensitivity / len(rna_names)
        # avg_ppv = sum_ppv / len(rna_names)
        # avg_mcc = sum_mcc / len(rna_names)

        # overall_sensitivity = 1.0 * sum_tp / (sum_tp + sum_fn)
        # overall_ppv = 1.0 * sum_tp / (sum_tp + sum_fp)
        tmp = 1.0 * (sum_tp + sum_fp) * (sum_tp + sum_fn) * (sum_tn + sum_fp) * (sum_tn + sum_fn)
        overall_MCC = (1.0 * sum_tp * sum_tn - 1.0 * sum_fp * sum_fn) / math.sqrt(tmp) * 100
        # overall_F1 = 2 * (overall_sensitivity * overall_ppv) / (overall_sensitivity + overall_ppv) * 100
        # overall_geometric_mean = math.sqrt(overall_sensitivity * overall_ppv) * 100

        f.write("ALL {:.1f}\n".format(overall_MCC))

        return overall_MCC


def sample_exp(n, lam):
    """Random sampling from an exponential distribution.

    Args:
        n (int): Number of data points to sample
        lam (float): Rate

    Returns:
        y (np.array): Sampled values

    """

    y = np.random.exponential(scale=1 / lam, size=n)

    return y


def sample_gev(n, mu, sigma, epsilon):
    """Random sampling from a generalized extreme value (GEV) distribution.

    Args:
        n (int): Number of data points to sample
        mu (float): Mean
        sigma (float): Variance
        epsilon (float): Shape

    Returns:
        y (np.array): Sampled values

    """

    epsilon *= -1  # Invert epsilon because scipy uses the inverse convention
    y = stats.genextreme.rvs(size=n, c=epsilon, loc=mu, scale=sigma)

    return y


def hmm_counter(structures, N=3):
    """Estimate HMM parameters based on RNA secondary structures encoded in numerical states.

    Args:
        structures (dict): Dictionary of structures encoded in numerical states.
        N (int): number of unique states.

    Returns:
        pi (np.array): Initial probabilities
        A (np.array): Transition probabilities
        rho (np.array): Ending probabilities

    """

    pi = np.zeros(N, dtype=float)
    # noinspection PyTypeChecker
    A = np.tile(0.0, (N, N))
    rho = np.zeros(N, dtype=float)
    n_cnt = 0

    for rna, structure in structures.items():
        structure = structure.astype(int)
        n = len(structure)

        # Initial
        pi[structure[0]] += 1

        # Transitions
        for i in range(n - 1):
            state_i = structure[i]
            state_j = structure[i + 1]
            A[state_i, state_j] += 1

        # Final
        rho[structure[-1]] += 1
        n_cnt += n

    # Normalize all probs
    pi /= np.sum(pi)
    A /= np.sum(A, axis=1)[:, np.newaxis]
    rho /= n_cnt

    return pi, A, rho


# noinspection PyTypeChecker
def hmm_generator(n, pi, A, rho, N=3):
    """Simulate paths based on HMM parameters.

    Args:
        n (int): Number of transcripts to simulate
        N (int): Number of unique states
        pi (np.array): Inital probabilities
        A (np.array): Transition probabilities
        rho (np.array): Ending probabilities
        N (int): Number of unique states

    Returns:
        transcripts (list): Simulated transcripts' paths

    """

    pool = np.arange(N)
    transcipts = []  # Hold all the simulated transcripts

    # Prepare rho - the ending probability
    rhop = np.tile(0.0, (N, 2))
    rho_pool = np.array([True, False])

    for i, p in enumerate(rho):
        rhop[i, 0] = p
        rhop[i, 1] = 1 - p

    for _ in range(n):

        path = []  # Current sequence path

        t = 0  # Current nucleotide position

        # t=0
        state = np.random.choice(pool, size=1, replace=True, p=pi)[0]  # Current state
        path.append(state)

        # t > 0
        while True:
            ended = np.random.choice(rho_pool, size=1, replace=True, p=rhop[state, :])[0]

            # Break if transcript stopped
            if ended:
                break

            state = np.random.choice(pool, size=1, replace=True, p=A[state, :])[0]
            path.append(state)

        transcipts.append(np.array(path, dtype=int))

    return transcipts


def path2shape(path, emit, N):
    """Path to SHAPE value using emission probability densities.

    Args:
        path (np.array): Sequence of states
        emit (dict): Emission probabilities parameters
        N (int): Number of unique states

    Returns:
        shape (np.array): Simulated shape

    """

    shape = np.repeat(np.nan, len(path))

    # Populate the emitted values for the current path
    for state in range(N):

        mask = path == state  # Select those states
        # noinspection PyTypeChecker
        n_draw = int(np.sum(mask))
        drew = np.nan

        if state == 0:  # Unpaired
            drew = sample_exp(n=n_draw,
                              lam=emit["unpaired"]["lambda"])
        elif state == 1:  # Helix-end
            drew = sample_gev(n=n_draw,
                              mu=emit["helix_end"]["mu"],
                              sigma=emit["helix_end"]["sigma"],
                              epsilon=emit["helix_end"]["epsilon"])
        elif state == 2:  # Stacked
            drew = sample_gev(n=n_draw,
                              mu=emit["stacked"]["mu"],
                              sigma=emit["stacked"]["sigma"],
                              epsilon=emit["stacked"]["epsilon"])
        # Replace the NaNs by the sampled values
        shape[mask] = drew

    return shape


def read_score(fp):
    scores = []
    with open(fp, "r") as f:
        next(f)  # skip header
        for line in f:
            line = line.rstrip().split()

            d = {"transcript": line[0],
                 "start": int(line[1]),
                 "end": int(line[2]),
                 "score": float(line[3]),
                 "cscore": float(line[4]),
                 "pvalue": float(line[4]),
                 "dot": line[5],
                 "path": [int(i) for i in list(line[6])],
                 "seq": line[7]}
            scores.append(d)
    return scores


def read_gammas(fp):
    rnas = {}
    with open(fp, "r") as f:

        while True:
            line1 = f.readline()
            line2 = f.readline()
            line3 = f.readline()

            if not line3:
                break

            rna = line1.split(">")[1].strip()
            rnas[rna] = np.array(line2.split(), dtype=float)
            # Line 3 is not used as it is 1 - line2

    return rnas


# noinspection PyUnboundLocalVariable
def read_n_shape(fp, n=np.inf):
    cnt = 0
    rnas = {}
    with open(fp, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                rna = line.split(">")[1].strip()
                cnt += 1
            else:
                rnas[rna] = np.array(line.split(), dtype=float)

            if cnt > n:
                break

    return rnas


def write_shape(fp, shape):
    """Write reactivities to a .shape formatted file."""

    with open(fp, "w") as f:
        for nuc in range(len(shape)):
            if np.isnan(shape[nuc]):
                f.write("{} {}\n".format(nuc + 1, -999))
            else:
                if shape[nuc] < 0:
                    f.write("{} {:.4f}\n".format(nuc + 1, 0))
                else:
                    f.write("{} {:.4f}\n".format(nuc + 1, shape[nuc]))


def write_fasta(fp, rna, seq):
    with open(fp, "w") as f:
        f.write("> {}\n {}".format(rna, seq))


def make_dir(path):
    """Create a directory. Spawn parents if needed."""

    if not os.path.exists(path):
        os.makedirs(path)


def print_attr(obj, attr=None):
    """Print all or selected attributes of an object to stdout. None prints all attributes."""

    for key in sorted(obj.__dict__):
        if attr is None:
            print(key.ljust(len(key) + 1, "\t"), " --> ", obj.__dict__[key])
        elif key in attr:
            print(key.ljust(len(key) + 1, "\t"), " --> ", obj.__dict__[key])


def rename_attribute(obj, old_name, new_name):
    """Rename a object attribute."""

    obj.__dict__[new_name] = obj.__dict__.pop(old_name)


def kwargs2attr(obj, kwargs):
    """Shallow copy attributes values based on keyword arguments. Assign only if the attribute already exists in obj."""
    if kwargs is None:
        return

    try:
        assert isinstance(kwargs, dict)
    except AssertionError:
        print("Attempted to set arguments not using a dictionary.")
        return

    attr = obj.__dict__.keys()
    for key in kwargs.keys():
        if key in attr:
            setattr(obj, key, kwargs[key])


def kwargs2attr_deep(obj, kwargs):
    """Deep copy attributes values based on keyword arguments. Assign only if the attribute already exists in obj."""
    if kwargs is None:
        return

    try:
        assert isinstance(kwargs, dict)
    except AssertionError:
        print("Attempted to set arguments not using a dictionary.")
        return

    attr = obj.__dict__.keys()
    for key in kwargs.keys():
        if key in attr:
            setattr(obj, key, deepcopy(kwargs[key]))


def file_length(fp):
    """Determine the number of rows in a file using linux wc -l."""

    cmd = subprocess.Popen(['wc', '-l', fp],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    result, err = cmd.communicate()
    if cmd.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def seconds_to_hms(t):
    """Formats a datetime.timedelta object to a HH:MM:SS string."""

    hours, remainder = divmod(t.total_seconds(), 3600)
    minutes, seconds = divmod(remainder, 60)
    hours, minutes, seconds = int(hours), int(minutes), int(seconds)
    if hours < 10:
        hours = "0%s" % int(hours)
    if minutes < 10:
        minutes = "0%s" % minutes
    if seconds < 10:
        seconds = "0%s" % seconds
    return "%s:%s:%s" % (hours, minutes, seconds)


def timer_start():
    """Start a timer."""

    return time.time()


def timer_stop(t0):
    """Stop a timer and return time as a formatted string."""

    t = timedelta(seconds=time.time() - t0)
    return seconds_to_hms(t)


def make_batches(n, batch_size, stochastic=False):
    """Creates a matrix for selecting batches.

    Args:
        n (int): Size of the dataset.
        batch_size (int): Size of each batch.
        stochastic (bool): Shuffle the data before making batch?

    Returns:
        sel (np.array): Array with each column being a boolean mask for selecting a single batch.

    """
    n = int(n)
    n_batches = int(np.ceil(n / batch_size))
    # noinspection PyTypeChecker
    sel = np.tile(False, (n, n_batches))

    # Make batch selectors
    ix = 0
    for i in range(n_batches):
        sel[ix:(ix + batch_size), i] = True
        ix += batch_size

    if stochastic:
        # Randomly shuffle the rows in the array
        np.random.shuffle(sel)

    return sel


def summarize_config(args):
    """Text summary of config arguments."""
    hline = "========================================================="
    text = "\n{}\n" \
           "Running pattern mining with the following parameters:\n" \
           "{}\n".format(hline, hline)

    for key in sorted(args.__dict__):
        text += "{}: {}\n".format(key, args.__dict__[key])

    text += hline

    return text


def selective_dict_deepcopy(input_dict, included_keys=None, exluded_keys=None):
    """Deep copy a dictionary using either keys to be included or keys to be excluded. If both are None, the entire
    dictionary is copied.

    Args:
        input_dict (dict): Input dictionary.
        included_keys (list): List of keys to include in the output dictionary.
        exluded_keys (list): List of keys to exclude in the output dictionary.

    Returns:
        Deep copy of the input dictionary with selected keys.

    """

    new_dict = {}
    keys_list = list(input_dict.keys())

    if included_keys is not None:
        keys_list = included_keys

    elif exluded_keys is not None:
        for exluded_key in list(exluded_keys):
            keys_list.remove(exluded_key)

    for key in keys_list:
        new_dict[key] = deepcopy(input_dict[key])

    return new_dict
