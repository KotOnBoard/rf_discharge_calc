# Set terminal and output file (e.g., PNG)
set terminal pngcairo size 1280,1024
set output 'multiplot_from_csv.png'

set datafile separator ","
set datafile columnheaders
set xdata time
set timefmt "%H:%M:%S" # Example: "YYYY-MM-DD HH:MM"
unset xlabel # For plots that are not at the bottomset xlabel "time [hh:mm]"
unset ylabel
set grid
set xtics rotate format ""

# Enable multiplot mode
# layout rows,columns arranges plots in a grid
set multiplot layout 4,1  spacing screen 1, 0.05 # Example: 2 rows, 2 columns

# Pressure
#set title "Plot 1: Data from Column 2"
#set ylabel "p [Pa]"
set bmargin 1 # No bottom margin for this plot (to align with the next plot's top)
plot 'testexport.csv' using "time":"p" with lines lc "red" lw 2 title "p [Pa]"

# Plot 2
set tmargin 0
#set yrange [60:80]
set offsets graph 0, 0, 0.2, 0.2
plot 'testexport.csv' using "time":"f" with lines lw 2 title "f0 [MHz]"

# Plot 3
#set yrange [0:*]
set offsets graph 0, 0, 2, 0.2
plot 'testexport.csv' using "time":"pr" with lines lw 2 title "Pr [W]",\
     'testexport.csv' using "time":"pf" with lines lw 2 title "Pf [W]"
     


# Plot 4
set xtics format "%H:%M"
set xlabel "time [hh:mm]"
set bmargin
plot 'testexport.csv' using "time":"urf" with lines lw 2 title "Urf [V]", \
     'testexport.csv' using "time":"ub" with lines lw 2 title "|Ub [V]|"

# Disable multiplot mode
unset multiplot