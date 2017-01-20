# Generate graphs comparing radial and division models with 1,2, or 3 parameters
# To generate the data, you must launch previously:
# . .../scripts/generate_data.sh radial_odd 3,5,7
# . .../scripts/generate_data.sh division_even 2,4,6

plot "radial_odd3rs.dat" using 1 title 'radial 1 parameter' with lines, \
 "radial_odd3rms.dat" using 1 lt 0 lc 1 notitle with lines, \
 "division_even2rs.dat" using 1 lc 2 title 'division 1 parameter' with lines, \
 "division_even2rms.dat" using 1 lt 0 lc 2 notitle with lines

plot "radial_odd5rs.dat" using 1 title 'radial 2 parameters' with lines, \
 "radial_odd5rms.dat" using 1 lt 0 lc 1 notitle with lines, \
 "division_even4rs.dat" using 1 lc 2 title 'division 2 parameters' with lines, \
 "division_even4rms.dat" using 1 lt 0 lc 2 notitle with lines

plot "radial_odd7rs.dat" using 1 title 'radial 3 parameters' with lines, \
 "radial_odd7rms.dat" using 1 lt 0 lc 1 notitle with lines, \
 "division_even6rs.dat" using 1 lc 2 title 'division 3 parameters' with lines, \
 "division_even6rms.dat" using 1 lt 0 lc 2 notitle with lines
