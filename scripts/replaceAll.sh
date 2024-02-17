#!/bin/bash
#*******************************************************************************
#****e* /replaceAll.sh
# NAME
# replaceAll.sh
# PURPOSE
# replace all occurences of X by Y in all files (using directory objects/)
#
# Please note, that this is a simple minded replacement and not a fully fledged
# refactoring tool. grep and sed will find pattern X, even if it is in a larger
# pattern.
#*******************************************************************************

if [ $# -ne 2 ]; then
    echo "please specify 2 command line arguments"
    exit 1
fi

PAT1=$1
PAT2=$2

echo ">${PAT1}<"
echo ">${PAT2}<"

FILES=`grep -l -i ${PAT1} objects/*90`

for file in $FILES; do
#    echo $file
    targ=`readlink $file`
    echo $targ
    sed -i "s|${PAT1}|${PAT2}|gi" "$targ"
done
