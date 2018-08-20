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
   3.2550000485032795        1340.0000000000000     
   6.5100000970065590        609.00000000000000     
   9.7650001455098394        285.00000000000000     
   13.020000194013118        116.00000000000000     
   16.275000242516398        51.000000000000000     
   19.530000291019679        22.000000000000000     
   22.785000339522959        16.000000000000000     
   26.040000388026236        7.0000000000000000     
   29.295000436529516        0.0000000000000000     
   32.550000485032797        0.0000000000000000     
e
   3.2550000485032795        1340.0000000000000     
   6.5100000970065590        609.00000000000000     
   9.7650001455098394        285.00000000000000     
   13.020000194013118        116.00000000000000     
   16.275000242516398        51.000000000000000     
   19.530000291019679        22.000000000000000     
   22.785000339522959        16.000000000000000     
   26.040000388026236        7.0000000000000000     
   29.295000436529516        0.0000000000000000     
   32.550000485032797        0.0000000000000000     
e
