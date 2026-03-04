# Set terminal and output file (e.g., PNG)
set terminal pngcairo size 1280,1024
set output './png/all_vs_ppl.png'

set datafile separator ","
#set datafile columnheaders
set grid xtics mxtics ytics mytics
set xtics 100
set mxtics 2
set key maxcols 2
set key maxrows 3
set termoption font "Times,18"

set multiplot title "p=3 [Pa] f=80 [MHz] S1=0.01 [m^2] l=0.1 [m]" \
	layout 2,2 spacing screen 1, 0.05 # Multiplot: 2 rows, 2 columns \
#set multiplot layout 2,2 spacing screen 1, 0.05 # Multiplot: 2 rows, 2 columns

unset xlabel
set ylabel "U_{bias} [V]"
set ytics 200
set mytics 4
set key right center
set object 1 rectangle from 275, graph 0 to 500, graph 1 fc rgb "red" fs solid 0.3 behind
set label 1 "U_{RF} > 500 [V]" at 400,0 center front font ",14"
plot './csv/all_vs_ppl.csv' using "ppl":"ub_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "ppl":"ub_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "ppl":"ub_he" with lines dt 3 lw 4 lc "magenta" title "He"

unset xlabel
set ylabel "U_{RF} [V]"
set key right center
set label 1 "U_{RF} > 500 [V]" at 400,200 center front font ",14"
plot './csv/all_vs_ppl.csv' using "ppl":"urf_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "ppl":"urf_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "ppl":"urf_he" with lines dt 3 lw 4 lc "magenta" title "He"

set key right bottom
set xlabel "Ppl [W]"
set ylabel "s [1]"
set yrange [-1:1]
set ytics 0.5
set mytics 4
set zeroaxis lc rgb "black" lw 2
set label 1 "U_{RF} > 500 [V]" at 400,0.75 center front font ",14"
plot './csv/all_vs_ppl.csv' using "ppl":"s_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "ppl":"s_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "ppl":"s_he" with lines dt 3 lw 4 lc "magenta" title "He"

set key right top
set ylabel "v_B [nm/min]"
set yrange [*:*]
set ytics auto
set label 1 "U_{RF} > 500 [V]" at 400,40 center front font ",14"
plot './csv/all_vs_ppl.csv' using "ppl":"d1h_b_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "ppl":"d1h_b_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "ppl":"d1h_b_he" with lines dt 3 lw 4 lc "magenta" title "He"

unset multiplot