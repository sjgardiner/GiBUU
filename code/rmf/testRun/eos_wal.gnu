set term postscript eps enhanced color "Helvetica" 40 lw 2
set out "eos_wal.eps"

set size 1.6,2.0
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

set xrange [0.:0.6]
set xlabel "{/Symbol r}_B (fm^{-3})"
set xtics 0.,0.1,0.6 scale 2
set mxtics 10

set yrange [-20.:130.]
set ylabel "E/A-m_N (MeV)"
set ytics -100.,20.,600. scale 2
set mytics 10

set key at 0.51,126. samplen 5 spacing 1.5 Right noreverse font "Helvetica,35"

plot "RMF_set10.dat" u 1:($5*1000.) t "Lang, NL1 (K=380 MeV, m*/m=0.83)" w l ls 1,\
     "RMF_set3.dat" u 1:($5*1000.) t "Lang, NL2 (K=210 MeV, m*/m=0.83)" w l ls 2,\
     "RMF_set11.dat" u 1:($5*1000.) t "Lang, NL3 (K=380 MeV, m*/m=0.70)" w l ls 3

set size
set autoscale
set ytics
unset label


