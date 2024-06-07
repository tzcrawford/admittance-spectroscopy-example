#!/bin/sh
maindir=$(pwd -P)/data
arrfile=$maindir/CT_out.tsv
logfile=$maindir/CTworkup.log


#/usr/bin/gnuplot -p <<EOF
$(which gnuplot) -p <<EOF

## GUI/Frontend and Exection Options
#Note, you may need to run gnuplot with the -p flag to keep gui window from immediately closing
#set term wxt                                ###GTK API terminal, allows for use of cursor with ruler etc.
#set term qt                                ###Qt API terminal
#set term x11                               ###Xterm, default for Xwindow systems
#set term png                               ###png terminal, Best for immediately outputting to a file
#set term png transparent truecolor background rgb '#000000' size 1600,1200 font "Liberation Mono:style=Bold,36" ###png terminal with transparent background but white on any opaque areas. Size 800x600 pixels.
set term png truecolor background rgb '#FFFFFF' size 1600,1200 font "Liberation Mono:style=Bold,36"
#set term svg enhanced size 500,500         ###svg terminal 5in by 5in

###Output Options
set output "arrhenius.png"                          ###sets output file and extension
#set output "a.jpg"
#set output "a.svg"

###Input Options
#set datafile separator ","                 #for plotting data from a comma separated variable (csv) file
set datafile separator "\t"                #for plotting data from a tab separated variable (tsv) file

###Legend options
#unset key                                  #disables legend 
set key top left Left reverse #enables legend and where
#set key box lc rgb '#000000'                #sets box around legend
#set key opaque                              #sets key on opaque background such that plots don't draw over it
set key tc rgb '#000000'                    #sets color of writing within key
set key font ",32"                         #legend font and size
#set key spacing 1.5                        #sets how far apart the lines in the key are vertically
set key samplen 1                          #sets how wide the key samples are

set border lw 4 lc rgb '#000000'                 #sets plot border a specific color



###Primary Graph Options
set title "Arrhenius plot" tc rgb '#000000'                           #title
set ylabel "ln(f/T_0^2) [ln(Hz/K^2)]" tc rgb '#000000'
set xlabel "1/T [1/K]" tc rgb '#000000'
#set xrange [4.5e-3:5.5e-3]                     #x axis range
set xrange [0:5.5e-3]                     #x axis range
#set yrange [2.2:3.6]                           #y axis range
set format y "%.1tE%T"                     #scientific notation
set format x "%.1tE%T"                     #scientific notation
#set xtics 5e-4                                #major x tick spacing
#set ytics 0.5                                 #major y ticks per major tick
set mxtics 10                                #minor x ticks per major tick
set mytics 5                               #minor y ticks per major tick
set grid                                    #enables grid on major ticks


###Functions
set samples 10000                           #enables the number of points plotted for functions (make larger if graph not smooth)

##Some custom linestyles
#lw is for line width (thickness)
#pt is for point type (which shape the marker is). I believe there are only 15 to choose from
#lc rgb '#~~~~~~' declares the line color
set style line 99 lw 3 pt 1 pointsize 3 lc rgb '#FFFFFF'      # white 
set style line 1 lw 3 pt 1 ps 3 lc rgb '#000000'      # black
set style line 2 lw 3 pt 2 ps 3 lc rgb '#CC0000'      # red 
set style line 3 lw 3 pt 3 ps 3 lc rgb '#000099'      # blue 
set style line 4 lw 3 pt 4 ps 3 lc rgb '#006600'      # green 
set style line 5 lw 3 pt 5 ps 3 lc rgb '#FF8000'      # orange 
set style line 6 lw 3 pt 6 ps 3 lc rgb '#FF007F'      # pink 
set style line 7 lw 3 pt 7 ps 3 lc rgb '#00FFFF'      # cyan 
set style line 8 lw 3 pt 8 ps 3 lc rgb '#606060'      # gray 
set style line 9 lw 3 pt 9 ps 3 lc rgb '#331900'      # brown 
set style line 10 lw 3 pt 10 ps 3 lc rgb '#FFFF00'     # yellow 
set style line 11 lw 3 pt 11 ps 3 lc rgb '#33FF33'     # light green 
set style line 12 lw 3 pt 12 ps 3 lc rgb '#FFB266'     # light orange 
set style line 13 lw 3 pt 13 ps 3 lc rgb '#FF99CC'     # light pink 
set style line 14 lw 3 pt 1 ps 3 lc rgb '#9999FF'     # light blue 
set style line 15 lw 3 pt 2 ps 3 lc rgb '#C0C0C0'     # light gray 

###Your Plots Here 
f(x)=`grep -P "^f\(x\) = " $logfile | sed "s/f(x) = //g"`

plot \
"`echo $arrfile`" using 5:7:6:8 title "100-250K" ls 1 w xyerrorbars, \
f(x) w l ls 1 lw 2 notitle

EOF
