#!/usr/bin/env python3
import numpy as np
import argparse
import os.path
import random
import warnings

from . import __version__, fastStructure
from .io import parse_bed, parse_str, write_output


def run_structure(G, K, output_prefix, tol, prior, cv):
    G = np.require(G, dtype=np.uint8, requirements='C')
    logfilename = '%s.%d.log' % (output_prefix, K)
    return fastStructure.infer_variational_parameters(G, K, logfilename, tol, prior, cv)


def main():
    prog = os.path.splitext(os.path.basename(__file__))[0]
    parser = argparse.ArgumentParser(prog=prog, description=
        "A user-friendly version of FastStructure, a fast algorithm for "
        "inferring population structure from large SNP genotype data.")

    parser.add_argument("-v", "--version", action="version", version="ezstructure " + __version__)
    parser.add_argument("-K", "--populations", metavar="<int>", type=int, required=True,
        help="Number of populations")
    parser.add_argument("-i", "--input", metavar="<file>", type=str, required=True,
        help="Path to input file")
    parser.add_argument("-o", "--output", metavar="<prefix>", type=str, required=True,
        help="Path and filename prefix for output files")
    parser.add_argument("-t", "--tol", metavar="<float>", type=float, default=1e-6,
        help="Convergence criterion; default: 1e-6")
    parser.add_argument("-p", "--prior", choices=("simple", "logistic"), default="simple",
        help="Choice of prior; default: simple")
    parser.add_argument("-c", "--cv", metavar="<int>", type=int, default=0,
        help="Number of test sets for cross-validation, <=0 implies no CV step; default: 0")
    parser.add_argument("-f", "--format", choices=("bed", "str"), type=str, default="bed",
        help="Format of input file; default: bed")
    parser.add_argument("-F", "--full", action="store_true",
        help="To output all variational parameters; optional")
    parser.add_argument("-s", "--seed", metavar="<int>", type=int,
        help="Manually specify seed for random number generator; optional")

    # parse command-line options
    args = parser.parse_args()
    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
    if args.tol <= 0:
        raise ValueError("A non-positive value was provided as convergence criterion")
    if args.populations <= 0:
        raise ValueError("A non-positive value was provided for the number of populations")

    # load data
    if args.format == 'bed':
        G = parse_bed(args.input)
    else:
        G = parse_str(args.input)

    # run the variational algorithm
    Q, P, other = run_structure(G, args.populations, args.output,
        args.tol, args.prior, args.cv)

    # write out inferred parameters
    write_output(Q, P, other, args.populations, args.output, args.full)


if __name__=="__main__":
    # ignore warnings with these expressions
    warnings.filterwarnings('ignore', '.*divide by zero.*',)
    warnings.filterwarnings('ignore', '.*invalid value.*',)
    main()
