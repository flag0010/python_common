## common.py
# Maintained by Lex Flagel
# Core utility functions for sampling, bioinformatics, and quick stats

import random
from collections import Counter, deque
from bisect import bisect_right
from itertools import combinations
from statistics import mean
import math

### --- Stats and math functions ---


def weighted_sampler(pop_dict, k=1):
    """Randomly sample keys from a dictionary using values as weights.
    Example:
    >>> m = {'a':3, 'b':2, 'c':5}
    >>> samps = weighted_sampler(m, k=20)
    """
    vals = list(pop_dict)
    weights = [pop_dict[v] for v in vals]
    return random.choices(vals, weights=weights, k=k)


def choose(n, k):
    """Binomial coefficient (n choose k). Requires Python 3.8+"""
    return math.comb(n, k)


def sampler(pop, size, replacement=False):
    """Sample from a population with or without replacement."""
    return random.choices(pop, k=size) if replacement else random.sample(pop, size)


def rank(x):
    """Return ranks of elements in list x. Ties are resolved by order of appearance."""
    return [i[0] for i in sorted(enumerate(x), key=lambda pair: pair[1])]


def order(x):
    """Return indices that would sort the list x."""
    return sorted(range(len(x)), key=lambda i: x[i])


def count_all(xlist, proportions=False):
    """Count all elements in a list. Optionally return proportions."""
    counts = Counter(xlist)
    if proportions:
        total = sum(counts.values())
        return {k: v / total for k, v in counts.items()}
    return counts


### --- Bioinformatics helpers ---


def revcom(s):
    """Return reverse complement of a DNA sequence (ACGT only)."""
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(trans)[::-1]


def get_fasta(file_name):
    """Read a FASTA file into a dict of {name: sequence}."""
    out = {}
    with open(file_name, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                curr_seq = line[1:]
                out[curr_seq] = []
            else:
                out[curr_seq].append(line)
    return {k: "".join(v) for k, v in out.items()}


def get_fasta_buffer(file_name):
    """Buffered FASTA reader. Yields (name, sequence) tuples."""
    with open(file_name) as file_iter:
        current_seq = []
        current_name = None
        for line in file_iter:
            line = line.strip()
            if line.startswith(">"):
                if current_seq and current_name is not None:
                    yield (current_name, "".join(current_seq))
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            yield (current_name, "".join(current_seq))


### --- File reading ---


def get_file(filename, splitchar="NA", buffered=False):
    """Read file into list of lists (or generator if buffered=True)."""
    opener = open(filename)
    if splitchar == "NA":
        return (
            (line.strip().split() for line in opener)
            if buffered
            else [line.strip().split() for line in opener]
        )
    else:
        return (
            (line.strip().split(splitchar) for line in opener)
            if buffered
            else [line.strip().split(splitchar) for line in opener]
        )


### --- Misc utilities ---


def sort_dict_by_val(d, reverse=False):
    """Return list of (key, value) tuples sorted by value (then key)."""
    return sorted(d.items(), key=lambda kv: (kv[1], kv[0]), reverse=reverse)


def pairwise(li):
    """Generate all pairwise combinations of elements in list."""
    return combinations(li, 2)


class rank_list:
    """Maintains a sorted list of top `maxlen` elements."""

    def __init__(self, maxlen=5):
        self.x = []
        self.maxlen = maxlen

    def add(self, i):
        b = bisect_right(self.x, i)
        self.x.insert(b, i)
        self.x = self.x[: self.maxlen]


def sliding_window(
    x,
    window_size,
    output_interval=1,
    flip_start=False,
    summary_stat=mean,
    start_output_frac=1.0,
):
    """simple sliding window implementation
    x = input values you want to slide windows over
    window_size = size of window
    output_interval = the number of strides you want the window to take before outputing a value, min=1, max=window_size, default=1
    flip_start = logical, do you wan to flip the first window and tack it onto the beginning. default = False
                 i.e. if x = [0,1,2,3,4,5,6,7,8] and window size is 3, then with flip_start=True
                 x becomes [2,1,0,0,1,2,3,4,5,6,7,8]. This fills the sliding window with a mirror image of the beginning values
                 and is a decent way (sometimes) to avoid outputing nothing while the window is filling
    summary_stat = any func that calcs a summary stat on the window, default=mean, but you could use median, sd, etc.
    start_outputing_frac = how many values in the window before you start outputing values, default=1.0 meaning wait
                 till window_size values are in before outputting anything.  Must be between 0-1
    """
    start_output_frac = float(start_output_frac)
    if flip_start:
        x = list(reversed(x[:window_size])) + x
    d = deque(maxlen=window_size)
    idx = 0
    for i in x:
        d.append(i)
        idx += 1
        if idx == output_interval and len(d) >= window_size * start_output_frac:
            yield summary_stat(d)
        if idx == output_interval:
            idx = 0


def flatten(lst):
    """Flatten a list of lists."""
    return [item for sublist in lst for item in sublist]


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]
