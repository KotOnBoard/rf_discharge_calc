# Set terminal and output file (e.g., PNG)
set terminal pngcairo size 1280,1024
set output './png/all_vs_f.png'

set datafile separator ","
#set datafile columnheaders
set grid xtics mxtics ytics mytics
set xtics 10
set mxtics 2
set key maxcols 2
set key maxrows 3
set termoption font "Times,18"

set multiplot title "p=3 [Pa] Ppl=100 [W] S1=0.04625 [m^2] l=0.1 [m]" \
	layout 2,2 spacing screen 1, 0.05 # Multiplot: 2 rows, 2 columns \
#set multiplot layout 2,2 spacing screen 1, 0.05 # Multiplot: 2 rows, 2 columns

unset xlabel
set ylabel "U_{bias} [V]"
set ytics 200
set mytics 4
set key right center
set object 1 rectangle from 30, graph 0 to 50, graph 1 fc rgb "red" fs solid 0.3 behind
set label 1 "U_{RF} > 500 [V]" at 40,-200 center front font ",14"
plot './csv/all_vs_f.csv' using "f":"ub_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "f":"ub_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "f":"ub_he" with lines dt 3 lw 4 lc "magenta" title "He"

unset xlabel
set ylabel "U_{RF} [V]"
set key right center
set label 1 "U_{RF} > 500 [V]" at 40,350 center front font ",14"
plot './csv/all_vs_f.csv' using "f":"urf_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "f":"urf_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "f":"urf_he" with lines dt 3 lw 4 lc "magenta" title "He"

set key right bottom
set xlabel "f [MHz]"
set ylabel "s [1]"
set yrange [-1:1]
set ytics 0.5
set mytics 4
set zeroaxis lc rgb "black" lw 2
set label 1 "U_{RF} > 500 [V]" at 40,-0.8 center front font ",14"
plot './csv/all_vs_f.csv' using "f":"s_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "f":"s_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "f":"s_he" with lines dt 3 lw 4 lc "magenta" title "He"

set key right top
set ylabel "v_B [nm/min]"
set yrange [*:*]
set ytics auto
set label 1 "U_{RF} > 500 [V]" at 40,15 center front font ",14"
plot './csv/all_vs_f.csv' using "f":"d1h_b_ar" with lines dt 1 lw 4 lc "blue" title "Ar", \
     ''			 using "f":"d1h_b_ne" with lines dt 2 lw 4 lc "red" title "Ne", \
     ''			 using "f":"d1h_b_he" with lines dt 3 lw 4 lc "magenta" title "He"

unset multiplot