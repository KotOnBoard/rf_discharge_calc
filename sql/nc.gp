# Read data from CSV file
#set datafile separator comma

# Set terminal and output
set terminal pngcairo enhanced font "Times,18" size 1200,800
set output './png/nc_vs_k.png'

# Set multiplot layout
set multiplot layout 1,3

# Extract min and max values for each gas
ar_min = 0.2791
ar_max = 0.60685

ne_min = 0.78408
ne_max = 1.04296

he_min = 3.02114
he_max = 4.66803

p = 1000
q = 100

nc_min(x) = 1/x

print "Argon:   Min = ", ar_min, " Max = ", ar_max
print "Neon:    Min = ", ne_min, " Max = ", ne_max  
print "Helium:  Min = ", he_min, " Max = ", he_max

# Common settings
set style fill transparent solid 0.3
set xrange [1:2.5]
set xtics 0.5
set mxtics 4
set grid xtics mxtics ytics mytics
set yrange [1:10000]
set logscale y
set xlabel "Inhomogeneity k [1]"
unset key

# Plot 1: Argon (Ar)
set title "Argon (Ar) Ratio"
set ylabel "N_C [1]"
d2m = ar_min
plot './csv/nc_out.csv' using 1:3 with lines lc rgb "red" lw 2 dt 1 title "Max Ar", \
     '' using 1:2 with lines lc rgb "red" lw 2 dt 1 title "Min Ar", \
     '' using 1:3:2 with filledcurves fc rgb "red" title "Ar Range"

# Plot 1: Neon (Ne)
set title "Neon (Ne) Ratio"
unset ylabel
d2m = ar_min
plot './csv/nc_out.csv' using 1:5 with lines lc rgb "red" lw 2 dt 1 title "Max Ne", \
     '' using 1:4 with lines lc rgb "red" lw 2 dt 1 title "Min Ne", \
     '' using 1:5:4 with filledcurves fc rgb "red" title "Ne Range"

# Plot 1: Helium (He)
set title "Helium (He) Ratio"
d2m = ar_min
plot './csv/nc_out.csv' using 1:7 with lines lc rgb "red" lw 2 dt 1 title "Max He", \
     '' using 1:6 with lines lc rgb "red" lw 2 dt 1 title "Min He", \
     '' using 1:7:6 with filledcurves fc rgb "red" title "He Range"

unset multiplot
