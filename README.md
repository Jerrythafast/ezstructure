# EasyStructure

## Introduction

EasyStructure is a cleant-up and repackaged version [fastStructure](https://github.com/rajanil/fastStructure). It runs on Python 3 and
can be straightforwardly installed using `pip install ezstructure`. It does not depend on any
non-Python libraries. Compared to the original *fastStructure* program, this version is about 50%
slower when using the `--prior=logistic` option (but not with the default value of `simple`).

*fastStructure* is a fast algorithm for inferring population structure from large SNP genotype data.
It is based on a variational Bayesian framework for posterior inference and was written in Python2.x.


## Citation

Anil Raj, Matthew Stephens, and Jonathan K. Pritchard. *fastSTRUCTURE: Variational Inference of
Population Structure in Large SNP Data Sets*, (Genetics) June 2014 197:573-589
[[Genetics](https://www.genetics.org/content/197/2/573.full),
[Biorxiv](https://biorxiv.org/content/early/2013/12/02/001073)]


## Installation

EasyStructure can be installed easily by using the `pip` package installer.
If you have `pip` installed, just run the following command to install the latest release:

    pip install -U ezstructure


## Using EasyStructure

EasyStructure is command-line software. The main command is `ezstructure`. You can view the
command-line help by running the command:

    ezstructure --help

EasyStructure, like its ancestor FastStructure, performs inference for the simplest, independent-loci,
admixture model, with two choices of priors that can be specified using the `--prior` flag.
Thus, unlike Structure, EasyStructure does not require the mainparams and extraparam files.
The inference algorithm used by FastStructure is fundamentally different from that of Structure and
requires the setting of far fewer options. All options can be passed via the flags listed above.


### Main options

The key options to pass to the scripts are the input file, the output file and the number of populations.
Assuming the input file is named `genotypes.bed` (with corresponding `genotypes.fam` and `genotypes.bim`),
the output file is named `genotypes_output` and the number of populations you would like is 3,
you can run the algorithm as follows:

    ezstructure -K 3 --input=genotypes --output=genotypes_output

This generates a `genotypes_output.3.log` file that tracks how the algorithm proceeds, and files
`genotypes_output.3.meanQ` and `genotypes_output.3.meanP` containing the posterior mean of
admixture proportions and allele frequencies, respectively. The orders of samples and
SNPs in the output files match those in the `.fam` file and `.bim` file, respectively.
Note that input file names need not include suffixes (e.g., .bed).


### Input data format

The current implementation can import data from [plink bed](https://zzz.bwh.harvard.edu/plink/binary.shtml)
format and the original Structure format. If the data are in plink format, ensure that
bed, bim and fam files for the dataset are all present in the same path.

While the original Structure program allowed for a more flexible input format, fastStructure expects a more
specific Structure-like input format. Specifically, rows in the data file correspond to samples, with two rows per sample
(note that only diploids are handled by this software), and columns correspond to SNPs. The first 6 columns
of the file will be ignored; these typically would include IDs, metadata, etc. This software only
handles bi-allelic loci. The two alleles at each locus can be encoded as desired; however, missing data
should be encoded as -9.


## Running on test data

A test simulated dataset is provided in `test/testdata.bed` in the source repository at
[GitHub](https://github.com/Jerrythafast/ezstructure) with genotypes sampled for
200 individuals at 500 SNP loci. The output files in `test/` were generated as follows:

    ezstructure -K 3 --input=test/testdata --output=testoutput_simple --full --seed=100
    ezstructure -K 3 --input=test/testdata --output=testoutput_logistic --full --seed=100 --prior=logistic

Executing the code with the provided test data should generate a log file identical to the ones in `test/`,
(except for the numbers in the `Iteration_Time (secs)` column) as a final check that the source code
has been downloaded and compiled correctly. The algorithm scales linearly with number of samples,
number of loci and value of K; the expected runtime for a new dataset can be computed from the runtime in the above log file.


## Choosing model complexity

In order to choose the appropriate number of model components that explain structure in the dataset,
we recommend running the algorithm for multiple choices of K. We have provided a utility tool,
`structure_choosek`, to parse through the output of these runs and provide a reasonable range of
values for the model complexity appropriate for this dataset.

Assuming the algorithm was run on the test dataset for choices of K ranging from 1 to 10, and
the output flag was --output=test/testoutput_simple, you can obtain the model complexity
by doing the following:

    ezstructure_choosek --input=test/testoutput_simple

The output would look like:

    Model complexity that maximizes marginal likelihood = 2
    Model components used to explain structure in data = 4


## Visualizing admixture proportions

In order to visualize the expected admixture proportions inferred by EasyStructure, we have
provided a simple tool to generate [Distruct](https://web.stanford.edu/group/rosenberglab/distruct.html)
plots using the mean of the variational posterior distribution over admixture proportions.
The samples in the plot will be grouped according to population labels inferred by EasyStructure.
However, if the user would like to group the samples according to some other categorical label
(e.g., geographic location), these labels can be provided as a separate file using the flag --popfile.
The order of labels in this file (one label per row) should match the order of samples in the input data files.

Assuming the algorithm was run on the test dataset for K=5, and
the output flag was --output=test/testoutput_simple, you can generate a Distruct plot
by doing the following:

    ezdistruct -K 5 --input=test/testoutput_simple --output=test/testoutput_simple_distruct.svg


## Python interface

As EasyStructure can be installed using pip, it is possible to use it as a dependency for other
packages. To use EasyStructure from within Python code, use the following example:

    from ezstructure.io import parse_bed, parse_str, write_output
    from ezstructure.structure import run_structure

    # Parse input file.
    G = parse_bed("example.bed")  # Or parse_str("example.str")

    # Set parameters.
    K = 3
    out_prefix = "example"
    tol = 1e-6
    prior = "simple"
    cv = 0

    # Run algorithm.
    Q, P, other = run_structure(G, K, out_prefix, tol, prior, cv)

    # Write output.
    write_output(Q, P, other, K, out_prefix, full=True)
