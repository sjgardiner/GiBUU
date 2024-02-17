#!/bin/bash

for i in `ls ~/GiBUU/buuinput/hades/*.acc`; do
    ./PrintHAFT.x $i
    echo
done
