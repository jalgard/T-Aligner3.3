import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np


# INPUT : - reference sequence (cryptogene)
#         - alignments in TAF format
#         - figure size (reference start:end)
def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--r', '--fa',
        action='store',
        help="Reference fasta file")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output table name")

    optionParser.add_argument('--d', '--dest',
        action='store',
        help="File DE analysis results [EdgeR format]")

    optionParser.add_argument('--s', '--sep',
        action='store',
        help="Separator [default: ,]")

    return optionParser

def ReadFasta(input_lines, split_space = False):

    fasta_data = []

    for line in input_lines:
        line = line.rstrip()
        if len(line) < 1:
            continue
        if line[0] == '>':
            if split_space == False:
                fasta_data.append([line[1:],''])
            else:
                fasta_data.append([line[1:].split(' ')[0],''])
        else:
            fasta_data[-1][1] += line.lstrip().rstrip().upper()

    for i, entry in enumerate(fasta_data):
        fasta_data[i].append(len(entry[1]))

    return fasta_data


if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])

    inputFastaData = None
    if runArgs.r is not None:
        if os.path.isfile(runArgs.r):
            inputFastaData = open(runArgs.r, 'r').readlines()
    else:
        sys.stderr.write('No reference file found or -r option missing!')
        exit()

    fastaEntries = ReadFasta(inputFastaData)
    inputFastaName = 0

    referenceSequence = fastaEntries[inputFastaName][1]
    # get cryptogene reference length
    ref_tless_length = 0

    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            ref_tless_length += 1

    # fill matrices with zero values
    pval_matrix = []
    fc_matrix = []
    for y in range(-20, 20):
        pval_matrix.append([0.0 for x in range(ref_tless_length)])
        fc_matrix.append([0.0 for x in range(ref_tless_length)])
    # read TAF
    separator = ','
    if runArgs.s is not None:
        if os.path.isfile(runArgs.s):
            separator = runArgs.s

    inputDEFile = ''
    if runArgs.d is not None:
        if os.path.isfile(runArgs.d):
            inputDEFile = runArgs.d
    else:
        sys.stderr.write('No DE file found or -d option missing!')
        exit()

    cy = 0
    cx = 0
    de_data = open(inputDEFile, 'r').readlines()[1:]
    for line in de_data:
        toks = line.split(separator)

        if 'N' not in toks[6] and float(toks[6]) < 1e-3 and (2 ** float(toks[2]) > 3 or 2 ** float(toks[2]) < 0.33):
            pval_matrix[cy][cx] = float(toks[6])
            fc_matrix[cy][cx] = 2 ** float(toks[2])
        else:
            pval_matrix[cy][cx] = 0
            fc_matrix[cy][cx] = 0
        cx += 1
        if cx == len(fc_matrix[0]):
            cx = 0
            cy += 1
            if cy == 40:
                cy = 0


    xi = []
    yi = []
    zi = []
    for y in range(len(pval_matrix)):
        for x in range(len(pval_matrix[0])):
            xi.append(10.0 * x)
            yi.append(410 - 10.0 * y)
            zi.append(pval_matrix[y][x])

    plt.figure(figsize=[ref_tless_length / 12, 7])
    plt.plot([-10.0, 10.0 * ref_tless_length + 10.0],[210,210], c='#a1d76a', lw=10, alpha=0.3)
    plt.scatter(xi,yi, c=zi, cmap='Oranges', s=5.0, lw=1.0)

    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)
    plt.savefig(output_pic, dpi=300)
