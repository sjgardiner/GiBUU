reset
set colors classic
unset bars

set term push
set terminal postscript eps enhanced color 'Helvetica' 24 size 5,4 dl 1
#set terminal pdf enhanced color font 'Helvetica,20' size 5,4 dl 1
set output 'CLAS.eps'

#set log y
#set mytics 10
#set format y '10^{%T}'

FF="($7/(2*pi*1e6))"

VERSION="r`svnversion -n ../..`"
sVERSION="{/=20 ".VERSION."}"

set lmargin 0
set bmargin 0
set tmargin 0
set rmargin 0

set multiplot layout 3,4 scale 1,1 offset 0.5,0

set yrange [0:39]
set xrange [1.3:1.95]

set mytics 5
set xtics 1.2,0.2,1.9
set mxtics 2

set format x ""

set format y
set ylabel "  "
set label 1 "{/=18 Q^2=0.425 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_04_045_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.425)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.425)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set format y ""
set ylabel "  "

set label 1 "{/=18 Q^2=0.475 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_045_05_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.475)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.475)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set label 1 "{/=18 Q^2=0.525 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_05_055_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.525)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.525)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set label 1 "{/=18 Q^2=0.575 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_055_06_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.575)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.575)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""


set format y
set ylabel "{/Symbol s}({/Symbol g}* p {/Symbol \256} p {/Symbol p^+ p^-}) [{/Symbol m}b]"
set label 1 "{/=18 Q^2=0.625 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_06_065_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.625)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.625)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set format y ""
set ylabel "  "
set label 1 "{/=18 Q^2=0.675 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_065_07_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.675)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.675)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set label 1 "{/=18 Q^2=0.725 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_07_075_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.725)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.725)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set label 1 "{/=18 Q^2=0.775 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_075_08_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.775)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.775)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""



set format x
set xlabel "  "
set format y
set ylabel "  "
set label 1 "{/=18 Q^2=0.825 GeV^2}" at graph 0.92, graph 0.92 right
set label 2 "{/=14 E_e = 2.039 GeV}" at graph 0.1, 0.75
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_08_085_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.825)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.825)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set xlabel "W [GeV]" offset graph 0.5, 0
set format y ""
set ylabel "  "
set label 1 "{/=18 Q^2=0.875 GeV^2}" at graph 0.92, graph 0.92 right
set label 2 ""
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_085_09_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.875)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.875)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set xlabel "  "
set label 1 "{/=18 Q^2=0.925 GeV^2}" at graph 0.92, graph 0.92 right
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_09_095_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.925)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.925)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

set xlabel "  "
set label 1 "{/=18 Q^2=0.975 GeV^2}" at graph 0.92, graph 0.92 right
set label 2 sVERSION at graph 0.2, 0.7 tc rgb "cyan"
plot "~/GiBUU/buuinput/electro_twoPi/Fedotov_data_Q2_0425_0975/cr_sect_Q2_095_010_int_stat_plus_sys_err.dat" u 1:2:3 w yerr lt 1 lw 3 pt 7 ps 1 lc 1 t "",\
     "fort.1300" u ($1):(($2==0.975)?(@FF):1/0) w l lt 1 lw 5 lc 3 t "",\
     "fort.1302" u ($1):(($2==0.975)?(@FF):1/0) w l lt 1 lw 2 lc 3 dt 2 t ""

unset multiplot

set term pop
set output

! fixbb CLAS.eps