# Set terminal and output file (e.g., PNG)
set terminal pngcairo size 1280,1024
set output './png/all_vs_f.png'

set datafile separator ","
#set datafile columnheaders
set grid
set key maxcols 2
set key maxrows 3
set title "p=3 [Pa] Ppl=100 [W] S1=0.04625 [m^2] l=0.1 [m] assy=8"
#set border 1 lc rgb "blue"
#set border 4 lc rgb "red"

set multiplot layout 2,2 spacing screen 1, 0.05 # Multiplot: 2 rows, 2 columns

set xlabel "f [MHz]"

set ylabel "Ub [V]" textcolor rgb "red"
plot './csv/all_vs_f.csv' using "f":"ub_ar" with lines dt 1 lw 2 lc "blue" title "Ub Ar", \
     ''			 using "f":"ub_ne" with lines dt 2 lw 2 lc "blue" title "Ub Ne", \
     ''			 using "f":"ub_he" with lines dt 3 lw 2 lc "blue" title "Ub He"

set ylabel "Urf [V]"
plot './csv/all_vs_f.csv' using "f":"urf_ar" with lines dt 1 lw 2 lc "red" title "Urf Ar", \
     ''			 using "f":"urf_ne" with lines dt 2 lw 2 lc "red" title "Urf Ne", \
     ''			 using "f":"urf_he" with lines dt 3 lw 2 lc "red" title "Urf He"

set ylabel "s [1]"
plot './csv/all_vs_f.csv' using "f":"s_ar" with lines dt 1 lw 2 lc "blue" title "s Ar", \
     ''			 using "f":"s_ne" with lines dt 2 lw 2 lc "blue" title "s Ne", \
     ''			 using "f":"s_he" with lines dt 3 lw 2 lc "blue" title "s He"

set ylabel "s [1]"
plot './csv/all_vs_f.csv' using "f":"d1h_b_ar" with lines dt 1 lw 2 lc "red" title "v_B Ar", \
     ''			 using "f":"d1h_b_ne" with lines dt 2 lw 2 lc "red" title "v_B Ne", \
     ''			 using "f":"d1h_b_he" with lines dt 3 lw 2 lc "red" title "v_B He"

unset multiplot