reset
set colors classic
unset bars

set term push
set terminal postscript eps enhanced color 'Helvetica' 28 size 5,6 dl 1
#set terminal pdf enhanced color font 'Helvetica,20' size 5,4 dl 1
set output 'MAID.Q0350.eps'

set lmargin 0
set bmargin 0
set tmargin 0
set rmargin 0

set multiplot layout 2,1 scale 1,1 offset 0.5,0


set log y
set mytics 10
set format y '10^{%T}'

set xrange [1.2:2.05]
set yrange [1:500]

set mxtics 5
#set mytics 5

FFnull="($6/(2*pi*1e6))"
FFplus="($8/(2*pi*1e6))"

VERSION="r`svnversion -n ../..`"
sVERSION="{/=24 ".VERSION."}"



set label 2 ""
set label 3 "{/=30 Q^2 = 0.35 GeV^2}" at graph 0.05, 0.2
set label 1 "{/=24 E_e = 7.2 GeV}" at graph 0.1, 0.1

set format x ""

set format y
set ylabel "  "
set label 4 "{/=34 {/Symbol g}* p {/Symbol \256} p {/Symbol p^0}}" at graph 0.5, graph 0.92 right

plot "~/GiBUU/buuinput/electro_onePi/MAID/gammaP_pi0P.Q2_0350.dat" u ($1/1000):($3+$4) w l lt 1 lw 4 lc -1 t "",\
     "~/GiBUU/buuinput/electro_onePi/MAID/gammaP_pi0P.Q2_0350.dat" u ($1/1000):($3+$4) every 5 w p pt 7 ps 2 lc -1  t "MAID",\
     "fort.130" u 1:($2==0.350?@FFnull:1/0) w l lt 1 lw 6 lc 1 t "GiBUU",\
     "fort.132" u 1:($2==0.350?@FFnull:1/0) w l lt 1 lw 2 lc 1 dt 2 t "GiBUU Res",\

#     "~/GiBUU/buuinput/electro_onePi/MAID/gammaP_pi0P.Q2_0450.dat" u ($1/1000):($3+$4) every 5 w p pt 6 ps 2 lc -1  t "MAID",\


set label 1 ""

set label 3 ""

set label 2 sVERSION at graph 0.1, 0.1 tc rgb "cyan"

set format x
set xlabel "W [GeV]"
set format y
set ylabel "{/Symbol s}({/Symbol g}* p {/Symbol \256} N {/Symbol p}) [{/Symbol m}b]" offset 0, graph 0.5

set label 4 "{/=34 {/Symbol g}* p {/Symbol \256} n {/Symbol p^+}}" at graph 0.5, graph 0.92 right

plot "~/GiBUU/buuinput/electro_onePi/MAID/gammaP_pi+N.Q2_0350.dat" u ($1/1000):($3+$4) w l lt 1 lw 4 lc -1 t "",\
     "~/GiBUU/buuinput/electro_onePi/MAID/gammaP_pi+N.Q2_0350.dat" u ($1/1000):($3+$4) every 5 w p pt 7 ps 2 lc -1  t "",\
     "fort.130" u 1:($2==0.350?@FFplus:1/0) w l lt 1 lw 6 lc 1 t "",\
     "fort.132" u 1:($2==0.350?@FFplus:1/0) w l lt 1 lw 2 lc 1 dt 2 t "",\





#"~/GiBUU/workingCode/code/init/ElectronGenerator/testRun/fort.721" u 1:($5*1000) w l lt 1 lw 5 lc -1 dt 2 t ""

unset multiplot

set term pop
set output

! fixbb MAID.Q0350.eps
