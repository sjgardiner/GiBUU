set term postscript eps enhanced color "Helvetica" 40 lw 2
set out "m_pm.eps"

set size 1.6,1.9
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

set xrange [0.:2.0]
set xlabel "{/Symbol r}_B (fm^{-3})"
set xtics 0.,0.2,2.0 scale 2
set mxtics 10

set yrange [200.:1600.]
set ylabel "m_{/Symbol \261} (MeV)"
set ytics 200.,200.,1600. scale 2
set mytics 10

set key at 1.9,1550. samplen 5 spacing 1.5 Right noreverse font "Helvetica,35"

plot "RMF_set31_sum.dat" u ($3+$4):($7*1000.) t "Zschiesche, Set P3" w l ls 2,\
     "RMF_set31_sum.dat" u ($3+$4):($8*1000.) notitle w l ls 2,\
     "RMF_set34_sum.dat" u ($3+$4):($7*1000.) t "Shin, Set 2" w l ls 1,\
     "RMF_set34_sum.dat" u ($3+$4):($8*1000.) notitle w l ls 1,\
     "RMF_set3.dat" u 1:($4*1000.) t "Lang, NL2" w l ls 3,\
     "RMF_set3.dat" u 1:($4*1000.-938.+1535.) notitle w l ls 3


#     "RMF_set33_sum.dat" u ($3+$4):($7*1000.) t "Shin, Set 1" w l ls 2,\
#     "RMF_set33_sum.dat" u ($3+$4):($8*1000.) notitle w l ls 2,\
#     "RMF_set32_sum.dat" u ($3+$4):($7*1000.) t "Zschiesche, Set P2, K=374 MeV" w l ls 1,\
#     "RMF_set32_sum.dat" u ($3+$4):($8*1000.) notitle w l ls 1,\
#     "RMF_set2.dat" u 1:($4*1000.) t "Lalazissis, NL3, K=272 MeV" w l ls 3     



set size
set autoscale
set ytics
unset label 

