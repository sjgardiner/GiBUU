reset
set colors classic
unset bars

set term push
set terminal postscript eps enhanced color 'Helvetica' 28 size 5,4 dl 2
#set terminal pdf enhanced color font 'Helvetica,20' size 5,4 dl 1
set output 'testPotNucleus.eps'

#set log y
#set mytics 10
#set format y '10^{%T}'

set xrange [0:12]
set mxtics 2
set xlabel "r [fm]"


set mytics 5
set ylabel "U [MeV]"



set key right bottom

set label 1 "{/=20 thick solid:    p = 0}" at graph 0.03,0.95
set label 2 "{/=20 thin dashed: p = p_F}" at graph 0.03,0.89

set label 3 "{/=34 Cu}" at graph 0.75,0.6
set label 4 "{/=18 EQS = 5}" at graph 0.75,0.52


plot "fort.111" u 1:($2*1e3) w l lt 1 lc 3 lw 6 t "n",\
     "fort.111" u 1:($3*1e3) w l lt 1 lc 3 lw 3 dt 2 t "",\
     "fort.111" u 1:($4*1e3) w l lt 1 lc 2 lw 6 t "p",\
     "fort.111" u 1:($5*1e3) w l lt 1 lc 2 lw 3 dt 2 t "",\
     "fort.111" u 1:($6*1e3) w l lt 1 lc 5 lw 3 t "{/Symbol D}^{- } ",\
     "fort.111" u 1:($7*1e3) w l lt 1 lc 5 lw 2 dt 2 t "",\
     "fort.111" u 1:($8*1e3) w l lt 1 lc 3 lw 3 t "{/Symbol D}^0 ",\
     "fort.111" u 1:($9*1e3) w l lt 1 lc 3 lw 2 dt 2 t "",\
     "fort.111" u 1:($10*1e3) w l lt 1 lc 2 lw 3 t "{/Symbol D}^+ ",\
     "fort.111" u 1:($11*1e3) w l lt 1 lc 2 lw 2 dt 2 t "",\
     "fort.111" u 1:($12*1e3) w l lt 1 lc 4 lw 3 t "{/Symbol D}^{++}",\
     "fort.111" u 1:($13*1e3) w l lt 1 lc 4 lw 2 dt 2 t "",\
     0 w l lt -1 lc 0 dt 2 lw 1 t ""


set term pop
set output
