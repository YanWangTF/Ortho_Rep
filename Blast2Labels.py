#! /usr/bin/python

#to transfer the e-value to log(e-value) for comparison
#usage: python Blast2Labels.py blastoutput.txt 1e-5 > LabelsfromBlast.txt

from sys import argv, stderr, stdout
from math import log10


FILE = open(argv[1],'r')

THRESH = float(argv[2])

for LINE in FILE:
    LINE=LINE.split()
    EV = float(LINE[10])
    if EV <= THRESH:
        stdout.write(LINE[0]+"\t"+LINE[1]+"\t")
        if EV == 0:
            LOG_EV = 200
        else:
            LOG_EV = -log10(EV)
        stdout.write(str(LOG_EV)+"\n")
