# Set terminal and output file (e.g., PNG)
set terminal pngcairo size 1280,1024
set output './png/ub_vs_ppl_vs_f.png'

set datafile separator ","
#set datafile columnheaders
set title "p=3 [Pa] S1=0.01 [m^2] l=0.1 [m] assy=8"

set xlabel "Ppl [W]"
set ylabel "f [MHz]"


set multiplot layout 2,3 # 2 rows, 1 column

set dgrid3d 100,100,2
set contour base
set view map
#unset xlabel
#unset ylabel
#unset zlabel
unset key

set title "p=3 [Pa] S1=0.01 [m^2] l=0.1 [m] assy=8\nS_a_r [arb.]"
set zrange [-1:1]
splot "./csv/all_vs_ppl_vs_f.csv" using "ppl":"f":"s_ar" with pm3d

set title "p=3 [Pa] S1=0.01 [m^2] l=0.1 [m] assy=8\nS_N_e [arb.]"
set zrange [-1:1]
splot "./csv/all_vs_ppl_vs_f.csv" using "ppl":"f":"s_ne" with pm3d

set title "p=3 [Pa] S1=0.01 [m^2] l=0.1 [m] assy=8\nS_H_e [arb.]"
set zrange [-1:1]
splot "./csv/all_vs_ppl_vs_f.csv" using "ppl":"f":"s_he" with pm3d

set zrange [*:*]
set title "p=3 [Pa] S1=0.01 [m^2] l=0.1 [m] assy=8\nd1H_B _A_r [arb.]"
splot "./csv/all_vs_ppl_vs_f.csv" using "ppl":"f":"d1h_b_ar" with pm3d,\
"" using 1:2:($4 > 500 ? 9 : 1/0) with points pt 7 ps 0.5 lc rgb "#000000" title "urf_ar > 500"
set title "p=3 [Pa] S1=0.01 [m^2] l=0.1 [m] assy=8\nd1H_B _N_e [arb.]"
splot "./csv/all_vs_ppl_vs_f.csv" using "ppl":"f":"d1h_b_ne" with pm3d
set title "p=3 [Pa] S1=0.01 [m^2] l=0.1 [m] assy=8\nd1H_B _H_e [arb.]"
splot "./csv/all_vs_ppl_vs_f.csv" using "ppl":"f":"d1h_b_he" with pm3d

unset multiplot