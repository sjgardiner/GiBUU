reset
set colors classic
unset bars

set term push
set terminal postscript eps enhanced color 'Helvetica' 28 size 5,4 dl 1
#set terminal pdf enhanced color font 'Helvetica,20' size 5,4 dl 1
set output 'branchingInt.eps'

set xlabel "W [GeV]"
set ylabel "branching ratio, 0 < M_{D_{13}(1520)} < W"

set mxtics 5
set mytics 2

set arrow 2 nohead from 1.524,0.0 to 1.524,1.0 ls 1 lc -1 lw 3 dt 4

#set log y
#set mytics 10
#set format y '10^{%T}'

plot "fort.1003" u 1:2 w l lt 2 lw 2 dt 2 t "{/Symbol p}N ",\
     "fort.1003" u 1:3 w l lt 4 lw 2 dt 3 t "{/Symbol pD}_0",\
     "fort.1003" u 1:4 w l lt 3 lw 2 dt 4 t "{/Symbol pD}_2",\
     "fort.1003" u 1:5 w l lt 1 lw 6 t "{/Symbol r}N "

set output 'branchingInt2.eps'

set xrange [1.2:1.6]
set yrange [0:0.4]
set mytics 5

set key left Left reverse

set arrow 2 nohead from 1.524,0.0 to 1.524,0.4 ls 1 lc -1 lw 3 dt 4

replot

set term pop
set output
