import numpy as np
import copy as cp
import matplotlib.pyplot as plt

"""bundlelib is a useful collection of function for working with motif bundles in the context
of RNA structure analysis."""

def read_nulls(fp, empty_repo, motif_map):

    null_score_repo = empty_repo

    with open(fp, 'r') as f:

        while f:

            read = f.readline().strip()

            if read == '':
                break

            path = motif_map[read]

            if path in null_score_repo.keys():
                data = np.array(f.readline().strip().split(), dtype='float64')
                null_score_repo[path] = data

            if not path:
                break

    return null_score_repo


def compile_bundle_dict(score_file, empty_transcript, motif_map, transcript_wise=1, score_column=4, path_column=6):

    bundle_dict = {}  # initialize dict

    with open(score_file, "r") as f:

        f.readline()  # Header line

        if not transcript_wise:
            # Use placeholder 'all' to store data
            bundle_dict['all'] = cp.deepcopy(empty_transcript)

        for line in f:

            row = line.strip().split()
            transcript = row[0]

            if not transcript_wise:

                transcript = 'all'

                path = row[path_column]  # Get path
                path = motif_map[path]  # Map to canonical name for dict structure

                if row[score_column] == 'nan':
                    continue  # Ignore nans

                score = float(row[score_column])

                if path not in bundle_dict[transcript].keys():
                    bundle_dict[transcript][path] = [score]
                else:
                    bundle_dict[transcript][path].append(score)
                continue


            if transcript not in bundle_dict.keys():
                # append transcript to data
                bundle_dict[transcript] = cp.deepcopy(empty_transcript)  # dict wrapper necessary to prevent overwriting

            path = row[path_column]
            c_score = float(row[score_column])

            start = int(row[1])
            motif = motif_map[path]

            # If we read a start position that is too long for us to store, we append the list with NaNs
            if (start+1) > len(bundle_dict[transcript][motif]):
                zeros_to_add = (start+1)-len(bundle_dict[transcript][motif])
                bundle_dict[transcript][motif] = cp.deepcopy(bundle_dict[transcript][motif])+[np.nan for i in range(zeros_to_add)]
                # The plus operator is appending lists

            if np.isnan(bundle_dict[transcript][motif][start]).any():
                bundle_dict[transcript][motif][start] = [c_score]  # Create list
            else:
                bundle_dict[transcript][motif][start].append(c_score)

    return bundle_dict


def read_posteriors(fp):

    ratios = []
    p0 = []
    p1 = []

    title_done = False
    p0_done = False

    with open(fp, "r") as f:

        for line in f:
            if line[0] == ">":
                title_done = True
                continue
            if title_done:
                p0.extend([float(x) for x in line.strip().split()[80:][:-80]])
                p0_done = True
                title_done = False
                continue
            if p0_done:
                p1.extend([float(x) for x in line.strip().split()[80:][:-80]])
                title_done = False
                p0_done = False
                continue

    # plt.hist(np.log(p0), np.linspace(-20, 0, 100), alpha=0.5)
    # plt.hist(p1, np.linspace(0, 1, 100), alpha=0.5)
    ratios = np.log(np.divide(p0, np.subtract(1, p0)))
    # ratios2 = np.log(np.divide(p1, p0))
    # ratios = np.array(ratios)
    # ratios = ratios[ratios<100]
    # plt.hist(ratios)
    # x = np.mean(ratios)
    x = np.mean(ratios)
    return x


def read_ensemble(fp):
    """Read in an ensemble of motifs to test from an ensemble input file."""

    bundle_map = {}
    bundle_map_r = {}
    empty_transcript = {}

    with open(fp, "r") as f:
        for line in f:
            if not line.split():
                continue
            motif = dot2states(line.split()[1], as_string=True)
            mutant = line.split()[0]
            if mutant not in bundle_map.keys():
                empty_transcript[mutant] = []
            bundle_map[motif] = mutant
            bundle_map_r[mutant] = motif

    return empty_transcript, bundle_map, bundle_map_r


def dot2states(dot, as_string=False):
    """Translate a dot-bracket string in a sequence of numerical states"""

    dot = dot.replace(".", "0")  # Unpaired
    dot = dot.replace("(", "1")  # Paired
    dot = dot.replace(")", "1")  # Paired
    dot = dot.replace(">", "1")  # Paired (ct2dot symbols)
    dot = dot.replace("<", "1")  # Paired (ct2dot symbols)
    dot = dot.replace("{", "1")  # Paired (ct2dot symbols)
    dot = dot.replace("}", "1")  # Paired (ct2dot symbols)

    if as_string:
        dotl = dot
    else:
        dotl = np.array(list(dot), dtype=np.int8)

    return dotl
