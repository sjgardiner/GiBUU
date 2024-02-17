set term postscript eps enhanced color "Helvetica" 45 lw 2
set out "XsRatio.eps"

set size 1.6,1.8
set origin 0.,0.

#set rmargin 12.

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

set xrange [0.:4.]
set xlabel "{/Symbol r}_B/{/Symbol r}_0"
set xtics 0.,0.5,4.0 scale 2
set mxtics 10
  
#set yrange [0.:1.]
#set ytics 0.,0.1,1. scale 2
#set mytics 10

set yrange [1.e-02:1.]
set log y
set format y "10^{%L}"
set ytics scale 2
set mytics 10
set ylabel "{/Symbol s}@^{med}_{/Symbol NN  \256 ND}/{/Symbol s}@^{vac}_{/Symbol NN  \256 ND}"


set title "Lalazissis, NL3 (K=271 MeV, m*/m=0.60)"

set key at 4.0,0.9 samplen 5 spacing 1.5 Right noreverse font "Helvetica,35"

rho0=0.168

f(x)=-alpha*x**beta
alpha=1.4
beta=3.0

fit [0.:3.] f(x) "XsectionRatio_2gev_RMF_set\ 2.dat" u ($1/rho0):(log($5)) via alpha,beta

plot "XsectionRatio_1gev_RMF_set\ 2.dat" u ($1/rho0):5 t "E_{lab}=1 GeV" w l ls 2,\
     "XsectionRatio_2gev_RMF_set\ 2.dat" u ($1/rho0):5 t "E_{lab}=2 GeV" w l ls 1,\
     exp(-1.3*x) t "exp(-1.3{/Symbol r}_B/{/Symbol r}_0)" w l ls 3,\
     exp(-1.7*x) t "exp(-1.7{/Symbol r}_B/{/Symbol r}_0)" w l ls 4


set size
set autoscale
set ytics
unset label 

