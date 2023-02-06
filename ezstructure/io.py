#!/usr/bin/env python3
import numpy as np
import os.path
import struct

# maps plink binary represntation of genotypes to an unsigned integer
# missing values are coded by `3`.
BED_GENOMAP = {'00':0, '01':1, '11':2, '10':3}


def parse_bed(prefix):
    new_prefix, ext = os.path.splitext(prefix)
    if ext.lower() == ".bed":
        # Strip extension.
        prefix = new_prefix

    # number of individuals
    with open(prefix+'.fam', 'rt') as handle:
        for i,line in enumerate(handle):
            pass
        Nindiv = i+1

    # Number of bytes to read in at a time
    Nbytes = Nindiv//4+(Nindiv%4>0)*1

    # number of SNPs
    with open(prefix+'.bim', 'rt') as handle:
        for i,line in enumerate(handle):
            pass
        Nsnp = i+1

    tobit = lambda x: ''.join([bin(i)[2:].zfill(8)[::-1] for i in struct.unpack('<%sB'%Nbytes, x)])
    genotype = np.zeros((Nindiv,Nsnp),dtype='uint8')

    with open(prefix+'.bed', 'rb') as handle:

        # check if the file is a valid plink bed file
        line = handle.read(1)
        checkA = bin(struct.unpack('<B', line)[0])[2:].zfill(8)[::-1]
        line = handle.read(1)
        checkB = bin(struct.unpack('<B', line)[0])[2:].zfill(8)[::-1]
        line = handle.read(1)
        checkC = bin(struct.unpack('<B', line)[0])[2:].zfill(8)[::-1]

        if checkA!="00110110" or checkB!="11011000":
            raise ValueError("This is not a valid bed file")

        # parse the bed file
        for l in range(Nsnp):
            line = handle.read(Nbytes)
            bytestr = tobit(line)
            for n in range(Nindiv):
                genotype[n,l] = BED_GENOMAP[bytestr[2*n:2*n+2]]

    return genotype


def parse_str(prefix):
    if not os.path.splitext(prefix)[1]:
        prefix += '.str'

    # read in data from file
    with open(prefix, 'rt') as handle:
        loci = np.array([line.strip().split()[6:] for line in handle])

    N = loci.shape[0]//2
    L = loci.shape[1]

    genotype = np.zeros((N,L),dtype=np.uint8)

    for l,locus in enumerate(loci.T):
        alleles = list(np.unique(locus))

        # assumes that missing alleles are coded as -9
        if '-9' in alleles:
            alleles.remove('-9')

        # selecting major and minor alleles
        major = alleles[0]
        try:
            minor = alleles[1]
        except IndexError:
            minor = ''

        # maps structure represntation of genotypes to an unsigned integer
        genomap = { minor+minor : 2,
                    minor+major : 1,
                    major+minor : 1,
                    major+major : 0
                  }

        # convert genotypes to fastStructure representation
        for n in range(N):
            try:
                geno = genomap[''.join(locus[2*n:2*(n+1)])]
                genotype[n,l] = geno
            except KeyError:
                # missing values are coded by `3`.
                genotype[n,l] = 3

    return genotype


def write_output(Q, P, other, K, prefix, full=False):
    """
    write the posterior means and variational parameters
    to separate output files.
    """
    with open('%s.%d.meanQ'%(prefix,K), 'wt') as handle:
        handle.write('\n'.join(['  '.join(['%.6f'%i for i in q]) for q in Q])+'\n')
    with open('%s.%d.meanP'%(prefix,K), 'wt') as handle:
        handle.write('\n'.join(['  '.join(['%.6f'%i for i in p]) for p in P])+'\n')
    if full:
        with open('%s.%d.varQ'%(prefix,K), 'wt') as handle:
            handle.write('\n'.join(['  '.join(['%.6f'%i for i in q]) for q in other['varQ']])+'\n')
        with open('%s.%d.varP'%(prefix,K), 'wt') as handle:
            handle.write('\n'.join(['  '.join(['%.6f'%i for i in np.hstack((pb,pg))])
                for pb,pg in zip(other['varPb'],other['varPg'])])+'\n')
