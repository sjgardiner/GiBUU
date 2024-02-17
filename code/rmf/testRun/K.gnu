set term postscript eps enhanced color "Helvetica" 40 lw 2
set out "K.eps"

set size 1.2,1.8
set origin 0.,0.

set for [i=1:7] linetype i dt i

set style line 1 lt 1 lw 5 lc rgb "black"
set style line 2 lt 2 lw 5 lc rgb "blue" 
set style line 3 lt 3 lw 5 lc rgb "brown"
set style line 4 lt 4 lw 5 lc rgb "red"
set style line 5 lt 5 lw 5 lc rgb "magenta"
set style line 6 lt 6 lw 5 lc rgb "purple"
set style line 7 lt 7 lw 5 lc rgb "green"
set style line 8 lt 8 lw 5 lc rgb "yellow"
set style line 9 lt 9 lw 5 lc rgb "violet"
set style line 10 lt 10 lw 5 lc rgb "dark-blue"
set style line 11 lt 11 lw 5 lc rgb "dark-violet"
set style line 12 lt 12 lw 5 lc rgb "dark-red"
set style line 13 lt 13 lw 5 lc rgb "dark-green"


file="RMF_set2.dat"


set xrange [0.145:0.150]
set xlabel "{/Symbol r} (fm^{-3})"
set xtics 0.13,0.005,0.17

set yrange [-0.1:0.1]
set ylabel "P (MeV/fm^3)"

x0=0.150  # rho0, fm^-3
K=300.    # incompressibility, MeV

P(x)=K/9.*(x-x0)

#fit [0.159:0.161]  P(x) file u 1:($8*1000.) via K,x0

fit [0.145:0.150]  P(x) file u 1:($6*1000.) via K,x0 

plot file u 1:($6*1000.) t "Set 2 (NL3)" w p,\
     P(x) t "fit"









set size
set autoscale
set ytics
unset label 

