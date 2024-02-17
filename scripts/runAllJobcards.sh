#!/bin/bash
#*******************************************************************************
#****e* /runAllJobcards.sh
# NAME
# runAllJobcards.sh
# PURPOSE
# This script generates a subdirectory in testRun/ and then runs all jobards
# it can find in the directory testRun/jobCards
#
# NOTES
# in order to parallelize the execution of the jobs, you may call this script
# as e.g.
#    Ntask=4 ./scripts/runAllJobcards.sh
# to run 4 jobs at the same time.
#
# cf.
# https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop
#*******************************************************************************

path=`pwd`
rundir=`date +"RunAll_%Y-%m-%d_%H-%M-%S"`
pathrundir="${path}/testRun/${rundir}"
echo "Generate: ${pathrundir}"
mkdir $pathrundir
cd $pathrundir

task()
{
    jobbase=`basename $1 .job`
    echo "Task: ${jobbase} ..."
    mkdir ${jobbase}
    cd ${jobbase}
    ln -s ${path}/testRun/jobCards/022_ExternalSource.inp 022_ExternalSource.inp
    /usr/bin/time ../../GiBUU.x < $1 &> ../${jobbase}.rep
    cd ..
    echo "Task: ${jobbase} finished!"
}

open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    printf '%.3d' $? >&3
    )&
}

if [ "x$Ntask" == "x" ]; then
    Ntask=1
fi
echo "Parallelize with Ntask=$Ntask"
open_sem $Ntask

for job in `ls ${path}/testRun/jobCards/*.job`; do
    run_with_lock task "$job"
done
