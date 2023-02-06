#!/usr/bin/env python3
import numpy as np
#import matplotlib.pyplot as plot  # Imported when running, as it is not a dependency of ezstructure.
import argparse
import colorsys
import os.path

from . import __version__

def plot_admixture(admixture, population_indices, population_labels, title):
    try:
        import matplotlib.pyplot as plot
    except ImportError:
        raise ValueError("Please install 'matplotlib' before running this script.")

    N,K = admixture.shape
    colors = [colorsys.hsv_to_rgb(h,0.9,0.7) for h in np.linspace(0,1,K+1)[:-1]]
    text_color = 'k'
    bg_color = 'w'
    fontsize = 12

    figure = plot.figure(figsize=(5,3))

    xmin = 0.13
    ymin = 0.2
    height = 0.6
    width = 0.74
    indiv_width = width/N
    subplot = figure.add_axes([xmin,ymin,width,height])
    for spine in subplot.spines.values():
        spine.set_linewidth(0.001)

    for k in range(K):
        if k:
            bottoms = admixture[:,:k].sum(1)
        else:
            bottoms = np.zeros((N,),dtype=float)

        lefts = np.arange(N)*indiv_width
        subplot.bar(lefts, admixture[:,k], width=indiv_width, bottom=bottoms, facecolor=colors[k], edgecolor=colors[k], linewidth=0.4)

        subplot.axis([0, N*indiv_width, 0, 1])
        subplot.tick_params(axis='both', top=False, right=False, left=False, bottom=False)
        xtick_labels = tuple(map(str,['']*N))
        subplot.set_xticklabels(xtick_labels)
        ytick_labels = tuple(map(str,['']*K))
        subplot.set_yticklabels(ytick_labels)

    position = subplot.get_position()
    title_position = (0.5, 0.9)
    figure.text(title_position[0], title_position[1], title, fontsize=fontsize, \
        color='k', horizontalalignment='center', verticalalignment='center')

    for p,popname in enumerate(population_labels):
        indices = np.where(population_indices==p)[0]
        if indices.size>0:
            vline_pos = (indices.max()+1)*indiv_width
            subplot.axvline(vline_pos, linestyle='-', linewidth=0.2, c='#888888')
            label_position = (xmin+(2*indices.min()+indices.size)*0.5*indiv_width, ymin-0.01)
            figure.text(label_position[0], label_position[1], popname, fontsize=6, color='k', \
                horizontalalignment='right', verticalalignment='top', rotation=70)

    return figure


def get_admixture_proportions(prefix, K, popfile=None):

    # load admixture proportions
    with open('%s.%d.meanQ'%(prefix,K), 'rt') as handle:
        admixture = np.array([line.strip().split() for line in handle]).astype('float')
    N,K = admixture.shape
    admixture = admixture/admixture.sum(1).reshape(N,1)

    # get population labels
    if popfile is not None:
        with open(popfile, 'rt') as handle:
            populations = [line.strip() for line in handle]
        population_labels = list(set(populations))

        # group populations by cluster similarity
        population_cluster = [np.mean(admixture[[i for i,p in enumerate(populations) if p==label],:],0).argmax() \
            for label in population_labels]
        population_labels = [population_labels[j] for j in np.argsort(population_cluster)]
        population_indices = np.array([population_labels.index(pop) for pop in populations])

        # re-order samples in admixture matrix
        order = np.argsort(population_indices)
        population_indices = population_indices[order]
        admixture = admixture[order,:]

    else:
        print("creating population labels based on inferred admixture proportions")
        population_labels = ['population %d'%i for i in range(1,K+1)]
        population_indices = np.argmax(admixture,1)

        # re-order samples in admixture matrix
        order = np.argsort(population_indices)
        population_indices = population_indices[order]
        admixture = admixture[order,:]
        order = np.arange(N)
        for k in range(K):
            order[population_indices==k] = order[population_indices==k][np.argsort(admixture[population_indices==k,:][:,k])[::-1]]
        admixture = admixture[order,:]

    return admixture, population_indices, population_labels


def main():
    prog = os.path.splitext(os.path.basename(__file__))[0]
    parser = argparse.ArgumentParser(prog=prog, description="Distruct plot tool for Structure")
    parser.add_argument("-v", "--version", action="version", version="ezstructure " + __version__)
    parser.add_argument("-K", "--populations", metavar="<int>", type=int, required=True,
        help="Number of populations")
    parser.add_argument("-i", "--input", metavar="<prefix>", type=str, required=True,
        help="Path to and prefix of input file names; same as output argument of Structure")
    parser.add_argument("-o", "--output", metavar="<file>", type=str, required=True,
        help="Output file to write")
    parser.add_argument("-p", "--popfile", metavar="<file>", type=str,
        help="File with known categorical labels; optional")
    parser.add_argument("-t", "--title", metavar="<figure title>", type=str,
        help="A title for the figure; optional")

    # parse command-line options
    args = parser.parse_args()
    title = args.title if args.title is not None else args.input

    # get the data to be plotted
    admixture, population_indices, population_labels = get_admixture_proportions(
        args.input, args.populations, args.popfile)

    # plot the data
    figure = plot_admixture(admixture, population_indices, population_labels, title)
    figure.savefig(args.output, dpi=300)


if __name__=="__main__":
    main()
