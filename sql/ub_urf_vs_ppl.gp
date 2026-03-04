# Set terminal and output file (e.g., PNG)
set terminal pngcairo size 1280,1024
set output './png/ub_urf_vs_pwr.png'

set datafile separator ","
#set datafile columnheaders
set grid
set key maxcols 2
set key maxrows 3
set title "p=3 [Pa] f=80 [MHz] S1=0.01 [m^2] l=0.1 [m] assy=8"
#set border 1 lc rgb "blue"
#set border 4 lc rgb "red"

set xlabel "Ppl [W]"
set ylabel "Ubias [V]" textcolor rgb "blue"
set y2label "Urf [V]" textcolor rgb "red"
set ytics nomirror
set y2tics nomirror
plot './csv/all_vs_ppl.csv' using "ppl":"ub_ar" with lines dt 1 lw 2 lc "blue" title "Ub Ar", \
     ''			 using "ppl":"ub_ne" with lines dt 2 lw 2 lc "blue" title "Ub Ne", \
     ''			 using "ppl":"ub_he" with lines dt 3 lw 2 lc "blue" title "Ub He", \
     ''			 using "ppl":"urf_ar" axes x1y2 with lines dt 1 lw 2 lc "red" title "Urf Ar", \
     ''			 using "ppl":"urf_ne" axes x1y2 with lines dt 2 lw 2 lc "red" title "Urf Ne", \
     ''			 using "ppl":"urf_he" axes x1y2 with lines dt 3 lw 2 lc "red" title "Urf He"