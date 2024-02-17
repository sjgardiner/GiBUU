reset
set macros
set colors classic


set log y
set xrange [0:1.55]
set mxtics 2
set yrange [0.8:40000]
set mytics 10

unset bars

# all:
#=====
fOut="Joos76.A.eps"
! sort -n -k 1 -k 2 A/fort.1300 > /tmp/aaa
! sort -n -k 1 -k 2 B/fort.1300 > /tmp/bbb
! sort -n -k 1 -k 2 C/fort.1300 > /tmp/ccc

# only RES:
#==========
#fOut="Joos76.B.eps"
#! sort -n -k 1 -k 2 A/fort.1302 > /tmp/aaa
#! sort -n -k 1 -k 2 B/fort.1302 > /tmp/bbb
#! sort -n -k 1 -k 2 C/fort.1302 > /tmp/ccc



FF="($7/(2*pi*1e6))"

VERSION="r`svnversion -n ..`"
sVERSION="{/=40 ".VERSION."}"


set ylabel "{/Symbol s}({/Symbol g}* p {/Symbol \256} p {/Symbol p^+ p^-}) [{/Symbol m}b]"
set xlabel "Q^2 [GeV^2]"

set terminal push
set terminal postscript eps enhanced color 'Helvetica' 26 size 5.0,9.0 dl 2
set output fOut

set label 1 "{/=24 E_e = 7.2 GeV}" at graph 0.07, 0.205
set label 2 sVERSION at graph 0.6, 0.95 tc rgb "cyan"


set key left bottom Left reverse

set format y '10^{%T}'

plot "~/GiBUU/buuinput/electro_twoPi/Joos1976.data" i 0 u (($2+$1)/2):($3*1000):($2-$1)/2:($4*1000) w xyerror lt 1 lw 3 pt 7 ps 3 lc 1 t "W=1.35 {/=16 ( x1000)}",\
     "" i 1 u (($2+$1)/2):($3*100):($2-$1)/2:($4*100) w xyerror lt 1 lw 3 pt 6 ps 3 lc 2 t "W=1.45 {/=16 ( x100)}",\
     "" i 2 u (($2+$1)/2):($3*50):($2-$1)/2:($4*50) w xyerror lt 1 lw 3 pt 5 ps 3 lc 3 t "W=1.55 {/=16 ( x50)}",\
     "" i 3 u (($2+$1)/2):($3*20):($2-$1)/2:($4*20) w xyerror lt 1 lw 3 pt 4 ps 3 lc 4 t "W=1.65 {/=16 ( x20)}",\
     "" i 4 u (($2+$1)/2):($3*10):($2-$1)/2:($4*10) w xyerror lt 1 lw 3 pt 13 ps 3 lc 5 t "W=1.75 {/=16 ( x10)}",\
     "" i 5 u (($2+$1)/2):($3*5):($2-$1)/2:($4*5) w xyerror lt 1 lw 3 pt 12 ps 3 lc 7 t "W=1.90 {/=16 ( x5)}",\
     "" i 6 u (($2+$1)/2):($3*2):($2-$1)/2:($4*2) w xyerror lt 1 lw 3 pt 11 ps 3 lc 8 t "W=2.10 {/=16 ( x2)}",\
     "" i 7 u (($2+$1)/2):($3*1):($2-$1)/2:($4*1) w xyerror lt 1 lw 3 pt 10 ps 3 lc 9 t "W=2.35 {/=16 ( x1)}",\
     "/tmp/aaa" u ($2):(($1==1.35)?(@FF*1000):1/0) w l lt 1 lw 5 lc 1 t "",\
     "/tmp/aaa" u ($2):(($1==1.45)?(@FF*100):1/0) w l lt 1 lw 5 lc 2 t "",\
     "/tmp/aaa" u ($2):(($1==1.55)?(@FF*50):1/0) w l lt 1 lw 5 lc 3 t "",\
     "/tmp/aaa" u ($2):(($1==1.65)?(@FF*20):1/0) w l lt 1 lw 5 lc 4 t "",\
     "/tmp/aaa" u ($2):(($1==1.75)?(@FF*10):1/0) w l lt 1 lw 5 lc 5 t "",\
     "/tmp/bbb" u ($2):(($1==1.90)?(@FF*5):1/0) w l lt 1 lw 5 lc 7 t "",\
     "/tmp/bbb" u ($2):(($1==2.10)?(@FF*2):1/0) w l lt 1 lw 5 lc 8 t "",\
     "/tmp/ccc" u ($2):(($1==2.35)?(@FF*1):1/0) w l lt 1 lw 5 lc 9 t ""


set output
set terminal pop

#! fixbb @fOut
