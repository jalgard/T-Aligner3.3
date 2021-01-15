'''

    Read TAF file and build a graph

'''

import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
import networkx as nx

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="Reads mapped in TAF file format")

    optionParser.add_argument('--r', '--ref',
        action='store',
        help="Reference in FASTA file")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output files prefix")

    return optionParser

class ESNode:

    def __init__(self, tlS, upN, ES):
        self.tlS = tlS
        self.upN = upN
        self.ES = ES
        self.RS = False

    def append(self, sequence):
        return sequence + self.upN + self.ES

    def set_refstate(self):
        self.RS = True

    def ref_state(self):
        return self.RS

    def __hash__(self):
        return hash(str(self.tlS) + "_" + str(self.upN) + "_" + str(self.ES))

    def __eq__(self, node):
        return self.tlS == node.tlS and self.ES == node.ES

class ESEdge:

    def __init__(self):
        self.rs = 0

    def inc(self):
        self.rs += 1

    def dec(self):
        self.rs -= 1


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

def ParseReadString(tlS, seq, G, ref_state=False):

    data = []

    S = 0
    for p in range(len(seq)):

        upN = seq[p]
        ES = 0

        while p+1 < len(seq) and seq[p+1] == 'T':
            p  += 1
            ES += 1

        data.append([upN, ES, S])
        S += 1

    for u,v in zip(data, data[1:]):
        n1 = ESNode(tlS+u[2], u[0], u[1])
        n2 = ESNode(tlS+v[2], v[0], v[1])

        if ref_state:
            n1.set_refstate()
            n2.set_refstate()

        if G.has_node(n1) == False:
            G.add_node(n1)

        if G.has_node(n2) == False:
            G.add_node(n2)

        if G.has_edge(n1, n2):
            G[n1][n2]['object'].inc()

        else:
            G.add_edge(n1, n2, object=ESEdge())

if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])
    inputTafData = None
    inputRefFasta = None

    if runArgs.t is not None:
        if os.path.isfile(runArgs.t):
            inputTafData = open(runArgs.t, 'r').readlines()

    else:
        sys.stderr.write('No TAF file provided! Please add --taf key or check if the file exists.')
        exit()

    if runArgs.r is not None:
        if os.path.isfile(runArgs.r):
            inputRefFasta = ReadFasta(open(runArgs.r, 'r').readlines())[0][1]

    else:
        sys.stderr.write('No reference provided! Please add --ref key or check if the file exists.')
        exit()

    G = nx.Graph()
    ParseReadString(0, inputRefFasta, G, True)

    for line in inputTafData:
        toks = line.rstrip().split('\t')
        ParseReadString(int(toks[2]), toks[5], G)

    pos = {}
    for n in G.nodes():
        pos[n] = (10 * int(n.tlS), 10 * int(n.ES))

    for e in G.edges():
        if G.get_edge_data(e[0], e[1])['object'].rs < 10:
            G.remove_edge(e[0], e[1])

    colors  = ['red' if u.ref_state() and v.ref_state() else 'blue' for u, v in G.edges()]
    weights = [np.sqrt(G[u][v]['object'].rs) for u,v in G.edges()]
    nx.draw(G, pos, node_size=1, width=weights, edge_color=colors)
    plt.show()
