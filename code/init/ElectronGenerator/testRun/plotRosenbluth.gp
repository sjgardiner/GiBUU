reset

set colors classic

set xrange [1.1:2.2]

#set yrange [1:]
#set log y

#set yrange [0:]

set multiplot layout 3,1

set xlabel "W [GeV]"

set xtics 1.1,0.1,2.2
set mxtics 5


set ylabel "{/Symbol s} [{/Symbol m}b]"
set yrange [0:*]
plot "fort.222" u 1:3 w lp pt 7 t "L",\
     "fort.222" u 1:4 w lp pt 7 t "T",\
     "fort.122" u 1:3 w l lt 1 t "L",\
     "fort.122" u 1:4 w l lt 2 t "T"

set ylabel "{/Symbol s}_L [{/Symbol m}b]"
set yrange [0:*]
plot "fort.222" u 1:3 w lp pt 7 lt 1 t "L",\
     "fort.122" u 1:3 w l lt -1 lw 3 t "tot",\
     "fort.122" u 1:16 w l lt 3 t "Res",\
     "fort.122" u 1:17 w l lt 4 t "1pi",\
     "fort.122" u 1:18 w l lt 5 t "2pi",\
     "fort.122" u 1:19 w l lt 2 t "DIS",\

set ylabel "{/Symbol s}_T [{/Symbol m}b]"
set yrange [0:*]
plot "fort.222" u 1:4 w lp pt 7 lt 1 t "T",\
     "fort.122" u 1:4 w l lt -1 lw 3 t "tot",\
     "fort.122" u 1:24 w l lt 3 t "Res",\
     "fort.122" u 1:25 w l lt 4 t "1pi",\
     "fort.122" u 1:26 w l lt 5 t "2pi",\
     "fort.122" u 1:27 w l lt 2 t "DIS"
