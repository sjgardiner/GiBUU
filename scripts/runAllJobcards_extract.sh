#!/bin/bash
for i in `ls *rep`;
do
#       echo $i
        j=`basename -s ".rep" "$i"`
        o=`grep -c "BUU simulation: finished" $i`
        t=`awk '{ if (/avgtext/) { print $3; } }' $i`
        m=`du -BM $j| awk '{ print $1 }'`

	[ -z "$t" ] && t="--running--"

        echo $o $t $m '#' $j
done
