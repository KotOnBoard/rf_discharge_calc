# Set terminal and output file (e.g., PNG)
set terminal pngcairo size 1280,1024
set output './png/d1h_b_s_vs_f.png'

set datafile separator ","
#set datafile columnheaders
set grid
set key maxcols 2
set key maxrows 3
set title "p=3 [Pa] Ppl=100 [W] S1=0.01 [m^2] l=0.1 [m] assy=8"
#set border 1 lc rgb "blue"
#set border 4 lc rgb "red"

set xlabel "f [MHz]"
set ylabel "s [1]" textcolor rgb "blue"
set y2label "d1H_B [nm/min]" textcolor rgb "red"
set ytics nomirror
set y2tics nomirror
plot './csv/all_vs_f.csv' using "f":"s_ar" with lines dt 1 lw 2 lc "blue" title "s_B Ar", \
     ''			 using "f":"s_ne" with lines dt 2 lw 2 lc "blue" title "s_B Ne", \
     ''			 using "f":"s_he" with lines dt 3 lw 2 lc "blue" title "s_B He", \
     ''			 using "f":"d1h_b_ar" axes x1y2 with lines dt 1 lw 2 lc "red" title "d1H_B Ar", \
     ''			 using "f":"d1h_b_ne" axes x1y2 with lines dt 2 lw 2 lc "red" title "d1H_B Ne", \
     ''			 using "f":"d1h_b_he" axes x1y2 with lines dt 3 lw 2 lc "red" title "d1H_B He", \
	 ''			 using "f":"s2_ar" with lines dt 1 lw 2 lc "cyan" title "s_M_o_O_3 Ar", \
     ''			 using "f":"s2_ne" with lines dt 2 lw 2 lc "cyan" title "s_M_o_O_3 Ne", \
     ''			 using "f":"s2_he" with lines dt 3 lw 2 lc "cyan" title "s_M_o_O_3 He", \
     ''			 using "f":"d1h_b_ar" axes x1y2 with lines dt 1 lw 2 lc "pink" title "d1H_MoO3 Ar", \
     ''			 using "f":"d1h_b_ne" axes x1y2 with lines dt 2 lw 2 lc "pink" title "d1H_MoO3 Ne", \
     ''			 using "f":"d1h_b_he" axes x1y2 with lines dt 3 lw 2 lc "pink" title "d1H_MoO3 He"