#!/usr/bin/env python3
import numpy as np
import argparse
import glob
import os.path

from . import __version__
from .vars import utils

def parse_logs(files):

    """
    parses through log files to extract marginal
    likelihood estimates from executing the
    variational inference algorithm on a dataset.

    Arguments:

        files : list
            list of .log file names

    """

    marginal_likelihood = []
    for file in files:
        with open(file, 'rt') as handle:
            for line in handle:
                if 'Marginal Likelihood' in line:
                    m = float(line.strip().split('=')[1])
                    marginal_likelihood.append(m)
                    break

    return marginal_likelihood


def parse_varQs(files):

    """
    parses through multiple .meanQ files to extract the mean
    admixture proportions estimated by executing the
    variational inference algorithm on a dataset. This is then used
    to identify the number of model components used to explain
    structure in the data, for each .meanQ file.

    Arguments:

        files : list
            list of .meanQ file names

    """

    bestKs = []
    for file in files:
        with open(file, 'rt') as handle:
            Q = np.array([list(map(float,line.strip().split())) for line in handle])
            Q = Q/utils.insum(Q,[1])

        N = Q.shape[0]
        C = np.cumsum(np.sort(Q.sum(0))[::-1])
        bestKs.append(np.sum(C<N-1)+1)

    return bestKs


def main():
    prog = os.path.splitext(os.path.basename(__file__))[0]
    parser = argparse.ArgumentParser(prog=prog, description="Choose model complexity for Structure")
    parser.add_argument("-v", "--version", action="version", version="ezstructure " + __version__)
    parser.add_argument("-i", "--input", metavar="<prefix>", type=str, required=True,
        help="Path and filname prefix for input files; will read all <prefix>*.log and <prefix>*.meanQ files")
    args = parser.parse_args()

    files = glob.glob('%s.*.log'%args.input)
    Ks = np.array([int(file.split('.')[-2]) for file in files])
    marginal_likelihoods = parse_logs(files)

    bestKs = parse_varQs(glob.iglob('%s.*.meanQ'%args.input))

    print("Model complexity that maximizes marginal likelihood = %d"%Ks[np.argmax(marginal_likelihoods)])
    print("Model components used to explain structure in data = %d"%np.argmax(np.bincount(bestKs)))


if __name__=="__main__":
    main()
