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
   2.7300000406801699        0.0000000000000000     
   5.4600000813603398        27.000000000000000     
   8.1900001220405105        121.00000000000000     
   10.920000162720680        227.00000000000000     
   13.650000203400850        213.00000000000000     
   16.380000244081021        155.00000000000000     
   19.110000284761192        62.000000000000000     
   21.840000325441359        22.000000000000000     
   24.570000366121530        10.000000000000000     
   27.300000406801701        2.0000000000000000     
e
   2.7300000406801699        0.0000000000000000     
   5.4600000813603398        27.000000000000000     
   8.1900001220405105        121.00000000000000     
   10.920000162720680        227.00000000000000     
   13.650000203400850        213.00000000000000     
   16.380000244081021        155.00000000000000     
   19.110000284761192        62.000000000000000     
   21.840000325441359        22.000000000000000     
   24.570000366121530        10.000000000000000     
   27.300000406801701        2.0000000000000000     
e
