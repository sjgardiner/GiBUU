 set dgrid3d          30 ,          50
 set xlabel 'density[fm^-3]'
 set ylabel 'momentum[GeV]'
 set zlabel 'potential[GeV]'
 set dgrid3d 10,30
 set contour base
set title "meson #101"
splot "fort.101" w l
