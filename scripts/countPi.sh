#!/bin/bash

# script to print the inclusive cross sections for the final states
# pi+X, pi+pi+X, pi+pi+pi+X from analyzing a "OutChannels*dat" file
#
# This script does not work with mawk, but expects gawk to be installed

for i in `ls OutChannels*dat`;
do
    Egamma=`awk '$2 ~ /E_gamma/ {print $4}' $i`
    OnePi=`awk 'BEGIN { sum = 0 } $0 ~ /pi/ {sum += $1} END {print sum}' $i`
    TwoPi=`awk 'BEGIN { sum = 0 } $0 ~ /pi[0+-]\s+pi/ {sum += $1} END {print sum}' $i`
    ThreePi=`awk 'BEGIN { sum = 0 } $0 ~ /pi[0+-]\s+pi[0+-]\s+pi/ {sum += $1} END {print sum}' $i`
    total=`awk '$1 ~ /total/ {print $4}' $i`

    echo "$Egamma ${total} ${OnePi} ${TwoPi} ${ThreePi}"
done
