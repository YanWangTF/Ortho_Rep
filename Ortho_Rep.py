#! /usr/bin/python

#to extract representatives of each taxa from ortholog group suggested by MCL
#usage: python Ortho_Rep.py ConcatenatedSeqs.fasta MCLoutput LabelsfromBlast.txt 1e-5 MaxClusterSize AllowedMinimumNumberTaxa

import networkx as nx
from sys import argv, stderr, stdout
from os import mkdir,chdir
from itertools import combinations

def SPECIES(NAME):
    return NAME.split("_")[0]

FASTA_FILE = open(argv[1],'r')
MCL_FILE = open(argv[2],'r')
LABEL_FILE = open(argv[3],'r')
EV_THRESH = float(argv[4])# Minimum EV to allow;
MAX_CLUSTER = int(argv[5])
MIN_SPECIES = int(argv[6]) # Smallest cluster, after resolving/removing paralogs, to report

FASTA = {}

NETWORK = nx.Graph()

i=0

for LINE in MCL_FILE:
    i += 1
    if (i % 1E2) == 0:
        stderr.write("\rReading MCL File, Line %.1fK\t" % (i/1E3))
    GENES = LINE.split()
    if len(GENES) > MAX_CLUSTER:
        continue # Purely for computational reasons; clusters too big will explode runtime and max out RAM
    if len(GENES) < MIN_SPECIES:
        continue # If cluster doesn't have enough members to satisfy the min species argument then go to next iteration
    TAXA = { SPECIES(x) for x in GENES }
    if len(TAXA) < MIN_SPECIES:
        continue # Pass this cluster if the initial # of taxa is already too low
    NETWORK.add_edges_from(combinations(GENES,2),EV=0) # If all initial criteria pass, make a cluster of nodes/edges

stderr.write("\nDONE\n")

stderr.write("Finished with MCL; there are %s clusters in the network, and %s nodes.\n" % 
             (nx.number_connected_components(NETWORK),NETWORK.number_of_nodes()))

i=0

for LINE in LABEL_FILE:
    i += 1
    if (i % 1E4) == 0:
        stderr.write("\rReading Labels File, Line %.2fM\t" % (i/1E6))
    LINE = LINE.split()
    GENE1 = LINE[0]
    GENE2 = LINE[1]
    EV = float(LINE[2])
    if EV >= EV_THRESH:
        if NETWORK.has_edge(GENE1,GENE2):
            NETWORK[GENE1][GENE2]["EV"] = EV # If EV meets minimum and an edge exists between genes, add EV as weight

stderr.write("DONE\n")

ORTHOSETS = [] # List of sets to print


OUTDIR = "Output_EV-%s_MinSp-%s" % (EV_THRESH,MIN_SPECIES)

for i,CLUSTER in enumerate(nx.connected_component_subgraphs(NETWORK,copy=False)):
    if (i % 1E2) == 0:
        stderr.write("\rProcessing Cluster %.1fK\t" % (i/1E3))
    RESOLVED = [] # List of taxa without paralogs or with paralogs resolved
    TAXA = { SPECIES(x) for x in CLUSTER }
    for PAIR in combinations(CLUSTER,2):
        if SPECIES(PAIR[0]) == SPECIES(PAIR[1]):
            CLUSTER.remove_edge(*PAIR) # Remove edges between genes from the same taxa
    for TAXON in TAXA:
        PARALOGS = sorted([ (x,NETWORK.degree(x,"EV")) for x in CLUSTER if TAXON in x ],key=lambda x: x[1],reverse=True)
        RESOLVED.append(PARALOGS[0][0]) # Put the best paralog in the resolved list
    ORTHOSETS.append(RESOLVED) # Add the set of resolved genes, one per taxa, to the orthosets list

stderr.write("DONE\n")

NETWORK = None

i=0

for LINE in FASTA_FILE:
    i += 1
    if (i % 1E4) == 0:
        stderr.write("\rReading Fasta File, Line %.2fM\t" % (i/1E6))
    LINE = LINE.strip()
    if LINE and (LINE[0] == ">"):
        HEADER = LINE[1:].split()[0]
        FASTA[HEADER] = ""
    else:
        FASTA[HEADER] += LINE

stderr.write("DONE\n")

stderr.write("Printing Ortholog Sets.\n")

mkdir(OUTDIR) # Make output directory

chdir(OUTDIR) # Switch to output directory

OS_COUNT = 0
G_COUNT = 0

for i,OS in enumerate(sorted(ORTHOSETS,key=len,reverse=True)):
    OS_COUNT += 1
    OS_OUT = open("Orthoset_%s.fasta" % (i+1),'w') # Open an output file
    for GENE in OS:
        G_COUNT += 1
        SEQ = FASTA[GENE] # Retrieve the sequence
        OS_OUT.write(">"+GENE+"\n") # Print a FASTA header
        OS_OUT.write(SEQ+"\n") # Print the sequence

stderr.write("Finished.  Printed %s ortholog sets with %s total genes.\n" % (OS_COUNT,G_COUNT))
