set title 'Fluid Velocity'
set xlabel 'cell # along x-dimension'
set ylabel 'cell # along y-dimension'

set terminal png
set output 'final_state.png'

#set terminal postscript eps size 12,10 enhanced color font 'Verdana,30' linewidth 2
#set output 'final_state.eps'

plot 'final_state.dat' using 1:2:5 with image
