reset
set colors classic
unset bars
set macros


fOut="Joos76_W.eps"

set term push
set terminal postscript eps enhanced color 'Helvetica' 28 size 5,4 dl 1
#set terminal pdf enhanced color font 'Helvetica,20' size 5,4 dl 1
set output fOut

#set log y
#set mytics 10
#set format y '10^{%T}'

set xrange [1.2:3.0]
set yrange [0:50]

set mxtics 5
set mytics 5

FF="($7/(2*pi*1e6))"

VERSION="r`svnversion -n ..`"
sVERSION="{/=24 ".VERSION."}"

set ylabel "{/Symbol s}({/Symbol g}* p {/Symbol \256} p {/Symbol p^+ p^-}) [{/Symbol m}b]"
set xlabel "W [GeV]"

set label 1 "{/=24 E_e = 7.2 GeV}" at graph 0.7, 0.75
set label 2 sVERSION at graph 0.03, 0.95 tc rgb "cyan"

W(Egamma)=sqrt(0.938**2+2*0.938*Egamma)

plot "~/GiBUU/buuinput/electro_twoPi/Joos1976.data" u (($5+$6)/2):($1==0.3?$3:1/0):($6-$5)/2:($4) w xyerror lt 1 lw 3 pt 6 ps 3 lc 2 t "Q^2=0.35 GeV^2",\
     "" u (($5+$6)/2):($1==0.4?$3:1/0):($6-$5)/2:($4) w xyerror lt 1 lw 3 pt 5 ps 3 lc 1 t "Q^2=0.45 GeV^2",\
     "fort.1300" u 1:($2==0.350?@FF:1/0) w l lt 1 lw 5 lc 2 t "",\
     "fort.1300" u 1:($2==0.450?@FF:1/0) w l lt 1 lw 5 lc 1 t "",\
     "fort.1302" u 1:($2==0.350?@FF:1/0) w l lt 1 lw 2 dt 2 lc 2 t "",\
     "fort.1302" u 1:($2==0.450?@FF:1/0) w l lt 1 lw 2 dt 2 lc 1 t "",\



set term pop
set output

#! fixbb @fOut
