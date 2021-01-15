# convert TAF file to table suitable for DE analysis with EdgeR or DESeq2

import sys, argparse, os

# determine cryptogene length from cryptogene fasta

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--r', '--fa',
        action='store',
        help="Reference fasta file")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output table name")

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with read alignments in TAF format")

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
    es_matrix = []
    es_productive_matrix = []
    for y in range(-20, 20):
        es_matrix.append([0 for x in range(ref_tless_length)])
        es_productive_matrix.append([0 for x in range(ref_tless_length)])
    # read TAF

    inputTafFile = ''
    if runArgs.t is not None:
        if os.path.isfile(runArgs.t):
            inputTafFile = runArgs.t
    else:
        sys.stderr.write('No TAF alignments file found or -t option missing!')
        exit()

    with open(inputTafFile, 'r') as ifile:
        for line in ifile:
            toks = line.split('\t')
            s = int(toks[2])
            es = [int(x) for x in toks[4][:-1].split(';')[1:-2]]

            # skip VERY long insertions
            skip = False
            for ess in es:
                if ess > 20:
                    skip = True
                    break
            if skip:
                continue

            # detect if editing drops to reference (cutoff 7)
            for p, ess in enumerate(es):
                if s + p >= ref_tless_length:
                    continue
                es_matrix[20-ess][s + p] += 1
                if p > 7:
                    if es[p-8:p] == [0, 0, 0, 0, 0, 0, 0, 0]:
                        es_productive_matrix[20-ess][s + p] += 1


    # print DE table
    output_DE = 'taf_to_DE_output.txt'
    if runArgs.o is not None:
        output_DE = str(runArgs.o)

    ofile = open(output_DE, 'w')
    for y in range(len(es_matrix)):
        for x in range(len(es_matrix[y])):
            DE_site = '{}:{}'.format(x, y)
            ofile.writelines('{}\t{}\n'.format(DE_site, es_matrix[y][x]))
    ofile.close()
