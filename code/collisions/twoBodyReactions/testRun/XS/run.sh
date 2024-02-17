#!/bin/bash
#*******************************************************************************
#****e* /run.sh
# NAME
# run.sh
# PURPOSE
# This script starts the php Development Server in order to offer the web
# interface of the cross section plotter locally to the browser.
#
# NOTES
# as prerequisites you need:
# * php
# * gnuplot
# * CrossSectionPlotter.x
#*******************************************************************************

# check the prerequisites:

if [ "x$0" != "x${BASH_SOURCE[0]}" ]; then
    echo "You should not source this script. STOP!"
    return 1;
fi


if ! command -v gnuplot &> /dev/null
then
    echo "ERROR: 'gnuplot' could not be found. STOP!"
    exit
fi

if ! command -v php &> /dev/null
then
    echo "ERROR: 'php' could not be found. STOP!"
    exit
fi

file="../CrossSectionPlotter.x"
if [[ ! -x "$file" ]]
then
    echo "File '$file' is not executable or found. STOP!"
    exit
fi

echo
echo "Please open '127.0.0.1:8000' in your favorite browser..."
echo
echo

# start the server:
php -S 127.0.0.1:8000
