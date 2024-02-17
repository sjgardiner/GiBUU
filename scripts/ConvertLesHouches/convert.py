#!/usr/bin/env python3
"""
Convert output files in LesHouches format to a 'FinalEvents.dat' file.

This program reads all GiBUU output files in LesHouches format and writes
the content to a file, which is similar to the 'FinalEvents.dat' file
produced when using the neutrino init.

The code is based on the library 'pylhef' by Roberto Vidal,
https://github.com/jrvidal/pylhef
"""

import glob
from pylhef import *
from KF2ID import *


f = open("FinalEvents.dat", 'w')

print("# 1:Run  2:Event  3:ID 4:Charge     5:perweight   6:position(1)   7:position(2)   8:position(3)   9:momentum(0)  10:momentum(1)  11:momentum(2)  12:momentum(3)     13:history  14:production_ID  15:enu", file=f)

files = glob.glob("EventOutput.Pert.*.lhe")

iRun = 0
for fLhe in sorted(files):
    data = read(fLhe)
    iRun=iRun+1

    iEv = 0
    for ev in data.events:
        iEv += 1
        for part in ev.particles:
            ID,IQ = KF_CODES[part.id]
            print("%7d%7d%7d%7d%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%11d%4d%14.6e"
                  % ( iRun, iEv, ID, IQ, ev.weight, 0.,0.,0., part.p[0], part.p[1], part.p[2], part.p[3], 0, 0, 0.), file=f )
