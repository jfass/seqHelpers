#!/usr/bin/python

"""
This script is intended to take fasta file (or stream)
and graph the cumulative coverage (given a genome size)
versus read length (longest first) given all sequences.
The purpose would be to choose a length cutoff - e.g.
for PBcR - such that all reads longer than the cutoff
provide sufficient genomic coverage (say, 60x) for
PBcR to work well. In addition to the graph, one
can provide a coverage cutoff and receive the read
length cutoff that provides the specified coverage.

Examples:
$ python est_genome_cov.py -g 2000000 reads.fasta
$ cat *.fa | python est_genome_cov.py -g 5000 -
$ python est_genome_cov.py -g 100000 -c 60 reads.fa
$ python est_genome_cov.py -g 5000 -f
"""

import argparse
import subprocess
import matplotlib.pyplot as plt

def getLengths(readsFileName):
    # open and iterate through reads file, gleaning lengths
    readsFile = open(readsFileName, 'r')
    lines = readsFile.readlines()
    lengthList = []
    for line in lines:
        lengthList.append(len(line))
    lengthList = sorted(lengthList, reverse=True)
    return lengthList

def cumulativeCovg(lengthList, genomeSize):
    # iterate through length list, creating fold covg list
    foldCovgs = []
    cumulativeCovg = 0
    for length in lengthList:
        cumulativeCovg = cumulativeCovg + float(length) / genomeSize
        foldCovgs.append(cumulativeCovg)
    return foldCovgs

def printGraph(lengths, foldCvgs):
    # paste lengths and coverages together, plot
    plt.plot(lengths, foldCvgs, marker='o', linestyle='-', color='b')
    plt.grid(True)
    plt.xlabel('read length (bp)')
    plt.ylabel('cumulative coverage (bp) of reads longer than x')
    plt.title('Genome Coverage vs Minimum Read Length Cutoff')
    plt.show()
    return()

def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('reads', help='reads in fasta format; no newlines in sequence')
    parser.add_argument('-g', '--genomeSize', help='size of genome, in base pairs', type=int)
    parser.add_argument('-c', '--coverageCutoff', help='coverage cutoff', type=float)
    parser.add_argument('-f', '--filter', action='store_true', help='output filtered reads to STDOUT')
    args = parser.parse_args()
    # calculate length and cumulative coverage
    lengths = getLengths(args.reads)
    coverages = cumulativeCovg(lengths, args.genomeSize)
    printGraph(lengths, coverages)

if __name__ == '__main__':
    main()
