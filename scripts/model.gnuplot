# copy/paste into gnuplot
# set terminal pdf color
#set autoscale
set yrange [0.00000000000000001:0.01]
set log y
set format y "%L"
unset label
set xtic auto
set ytic auto
set xlabel "lens index"
set ylabel "average/max residual error (in log_10)"

set term pdf
set output "radial.pdf"
set key left top
plot "radial3s.dat" using 1 title 'order 3' with lines, \
     "radial5s.dat" using 1 title 'order 5' with lines, \
     "radial7s.dat" using 1 title 'order 7' with lines, \
     "radial9s.dat" using 1 title 'order 9' with lines, \
     "radial3ms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "radial5ms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "radial7ms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "radial9ms.dat" using 1 lt 0 lc 4 notitle with lines

set output "radialr.pdf"
set key right bottom
plot "radial3rs.dat" using 1 title 'order 3' with lines, \
     "radial5rs.dat" using 1 title 'order 5' with lines, \
     "radial7rs.dat" using 1 title 'order 7' with lines, \
     "radial9rs.dat" using 1 title 'order 9' with lines, \
     "radial3rms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "radial5rms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "radial7rms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "radial9rms.dat" using 1 lt 0 lc 4 notitle with lines


set output "poly.pdf"
set key left top
plot "poly3s.dat" using 1 title 'order 3' with lines, \
     "poly5s.dat" using 1 title 'order 5' with lines, \
     "poly7s.dat" using 1 title 'order 7' with lines, \
     "poly9s.dat" using 1 title 'order 9' with lines, \
     "poly3ms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "poly5ms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "poly7ms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "poly9ms.dat" using 1 lt 0 lc 4 notitle with lines

set output "polyr.pdf"
set key right bottom
plot "poly3rs.dat" using 1 title 'order 3' with lines, \
     "poly5rs.dat" using 1 title 'order 5' with lines, \
     "poly7rs.dat" using 1 title 'order 7' with lines, \
     "poly9rs.dat" using 1 title 'order 9' with lines, \
     "poly3rms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "poly5rms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "poly7rms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "poly9rms.dat" using 1 lt 0 lc 4 notitle with lines


set output "division.pdf"
set key left top
plot "division3s.dat" using 1 title 'order 3' with lines, \
     "division5s.dat" using 1 title 'order 5' with lines, \
     "division7s.dat" using 1 title 'order 7' with lines, \
     "division9s.dat" using 1 title 'order 9' with lines, \
     "division3ms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "division5ms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "division7ms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "division9ms.dat" using 1 lt 0 lc 4 notitle with lines

set output "divisionr.pdf"
set key right bottom
plot "division3rs.dat" using 1 title 'order 3' with lines, \
     "division5rs.dat" using 1 title 'order 5' with lines, \
     "division7rs.dat" using 1 title 'order 7' with lines, \
     "division9rs.dat" using 1 title 'order 9' with lines, \
     "division3rms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "division5rms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "division7rms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "division9rms.dat" using 1 lt 0 lc 4 notitle with lines


set output "rational.pdf"
set key left top
plot "rational3s.dat" using 1 title 'order 3' with lines, \
     "rational5s.dat" using 1 title 'order 5' with lines, \
     "rational7s.dat" using 1 title 'order 7' with lines, \
     "rational9s.dat" using 1 title 'order 9' with lines, \
     "rational3ms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "rational5ms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "rational7ms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "rational9ms.dat" using 1 lt 0 lc 4 notitle with lines

set output "rationalr.pdf"
set key right bottom
plot "rational3rs.dat" using 1 title 'order 3' with lines, \
     "rational5rs.dat" using 1 title 'order 5' with lines, \
     "rational7rs.dat" using 1 title 'order 7' with lines, \
     "rational9rs.dat" using 1 title 'order 9' with lines, \
     "rational3rms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "rational5rms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "rational7rms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "rational9rms.dat" using 1 lt 0 lc 4 notitle with lines

set output "radial_odd.pdf"
set key left top
plot "radial_odd3s.dat" using 1 title 'order 3' with lines, \
     "radial_odd5s.dat" using 1 title 'order 5' with lines, \
     "radial_odd7s.dat" using 1 title 'order 7' with lines, \
     "radial_odd9s.dat" using 1 title 'order 9' with lines, \
     "radial_odd3ms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "radial_odd5ms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "radial_odd7ms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "radial_odd9ms.dat" using 1 lt 0 lc 4 notitle with lines


set output "radial_oddr.pdf"
set key right bottom
plot "radial_odd3rs.dat" using 1 title 'order 3' with lines, \
     "radial_odd5rs.dat" using 1 title 'order 5' with lines, \
     "radial_odd7rs.dat" using 1 title 'order 7' with lines, \
     "radial_odd9rs.dat" using 1 title 'order 9' with lines, \
     "radial_odd3rms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "radial_odd5rms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "radial_odd7rms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "radial_odd9rms.dat" using 1 lt 0 lc 4 notitle with lines

set output "FOV.pdf"
set key right bottom
plot "FOV1s.dat" using 1 title 'order 1' with lines, \
     "FOV5s.dat" using 1 title 'order 5' with lines, \
     "FOV7s.dat" using 1 title 'order 7' with lines, \
     "FOV9s.dat" using 1 title 'order 9' with lines, \
     "FOV1ms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "FOV5ms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "FOV7ms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "FOV9ms.dat" using 1 lt 0 lc 4 notitle with lines

set output "FOVr.pdf"
set key right bottom
plot "FOV1rs.dat" using 1 title 'order 1' with lines, \
     "FOV5rs.dat" using 1 title 'order 5' with lines, \
     "FOV7rs.dat" using 1 title 'order 7' with lines, \
     "FOV9rs.dat" using 1 title 'order 9' with lines, \
     "FOV1rms.dat" using 1 lt 0 lc 1 notitle with lines, \
     "FOV5rms.dat" using 1 lt 0 lc 2 notitle with lines, \
     "FOV7rms.dat" using 1 lt 0 lc 3 notitle with lines, \
     "FOV9rms.dat" using 1 lt 0 lc 4 notitle with lines

set term wxt
