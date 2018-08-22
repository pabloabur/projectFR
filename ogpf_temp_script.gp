# ogpf libray
# Rev. 0.22 of March 9th, 2018
# Licence: MIT

# gnuplot global setting
set term wxt size 640,480 enhanced font "verdana,10" title "ogpf libray: Rev. 0.22 of March 9th, 2018"

# ogpf extra configuration
# -------------------------------------------
# color definitions
set style line 1 lc rgb "#800000" lt 1 lw 2
set style line 2 lc rgb "#ff0000" lt 1 lw 2
set style line 3 lc rgb "#ff4500" lt 1 lw 2
set style line 4 lc rgb "#ffa500" lt 1 lw 2
set style line 5 lc rgb "#006400" lt 1 lw 2
set style line 6 lc rgb "#0000ff" lt 1 lw 2
set style line 7 lc rgb "#9400d3" lt 1 lw 2

# Axes
set border linewidth 1.15
set tics nomirror

# grid
# Add light grid to plot
set style line 102 lc rgb "#d6d7d9" lt 0 lw 1
set grid back ls 102

# plot style
set style data linespoints

# -------------------------------------------

 
# plot scale
 
# Annotation: title and labels
set title "ISI histogram"
set xlabel "ISI"
set ylabel "count"
 
# axes setting

plot "-" notitle with impulses lw 2.5 \
, "-" notitle with points pt 7
   3.1050000462681053        1.0000000000000000     
   6.2100000925362107        26.000000000000000     
   9.3150001388043169        188.00000000000000     
   12.420000185072421        280.00000000000000     
   15.525000231340528        199.00000000000000     
   18.630000277608634        97.000000000000000     
   21.735000323876740        27.000000000000000     
   24.840000370144843        6.0000000000000000     
   27.945000416412949        4.0000000000000000     
   31.050000462681055        2.0000000000000000     
e
   3.1050000462681053        1.0000000000000000     
   6.2100000925362107        26.000000000000000     
   9.3150001388043169        188.00000000000000     
   12.420000185072421        280.00000000000000     
   15.525000231340528        199.00000000000000     
   18.630000277608634        97.000000000000000     
   21.735000323876740        27.000000000000000     
   24.840000370144843        6.0000000000000000     
   27.945000416412949        4.0000000000000000     
   31.050000462681055        2.0000000000000000     
e
