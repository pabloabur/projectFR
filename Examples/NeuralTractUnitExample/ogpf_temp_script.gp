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
   2.6050000388175247        0.0000000000000000     
   5.2100000776350495        8.0000000000000000     
   7.8150001164525751        100.00000000000000     
   10.420000155270099        208.00000000000000     
   13.025000194087625        226.00000000000000     
   15.630000232905150        145.00000000000000     
   18.235000271722676        84.000000000000000     
   20.840000310540198        38.000000000000000     
   23.445000349357723        17.000000000000000     
   26.050000388175249        5.0000000000000000     
e
   2.6050000388175247        0.0000000000000000     
   5.2100000776350495        8.0000000000000000     
   7.8150001164525751        100.00000000000000     
   10.420000155270099        208.00000000000000     
   13.025000194087625        226.00000000000000     
   15.630000232905150        145.00000000000000     
   18.235000271722676        84.000000000000000     
   20.840000310540198        38.000000000000000     
   23.445000349357723        17.000000000000000     
   26.050000388175249        5.0000000000000000     
e
