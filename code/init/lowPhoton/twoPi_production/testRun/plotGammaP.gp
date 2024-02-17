reset
set colors classic
unset bars

#set term push
#set terminal postscript eps enhanced color 'Helvetica' 28 size 5,4 dl 1
#set terminal pdf enhanced color font 'Helvetica,20' size 5,4 dl 1
#set output 'dummy.eps'

set multiplot layout 2,1

set log y
set mytics 10
set format y '10^{%T}'
set yrange [1:]


set xlabel "W [GeV]"
set ylabel "sigma [mub]"

set key right bottom

plot "gammaProton.dat" u 2:4 w l t "p pi+ pi-",\
     "gammaProton.dat" u 2:5 w l t "n pi+ pi0",\
     "gammaProton.dat" u 2:6 w l t "p pi0 pi0",\

unset log y
set yrange [0:*]
set format y

set ylabel "relative strength"

plot "gammaProton.dat" u 2:($4/($4+$5+$6)) w l t "p pi+ pi-",\
     "gammaProton.dat" u 2:($5/($4+$5+$6)) w l t "n pi+ pi0",\
     "gammaProton.dat" u 2:($6/($4+$5+$6)) w l t "p pi0 pi0"


unset multiplot

#set term pop
#set output
