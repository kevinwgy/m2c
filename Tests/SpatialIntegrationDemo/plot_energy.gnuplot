#=========================================================
# Energy Variation History
#
# Plots:
#   Zone i: m_i(t) - m_i(0)
#   Total : sum_i m_i(t) - sum_i m_i(0)
#
# Assumptions:
#   - Files are in results/
#   - Each file has 5 columns
#   - Time is column 2
#   - Energy is column 3
#   - First 2 lines are headers
#   - All files have identical time histories
#=========================================================

reset

#---------------------------------------------------------
# Output
#---------------------------------------------------------
set terminal pngcairo enhanced size 1200,800 font "Arial,16"
set output "totalenergy_variation.png"

#---------------------------------------------------------
# Appearance
#---------------------------------------------------------
set border linewidth 1.5

set tics out
set mxtics
set mytics

set grid xtics ytics mxtics mytics \
    lw 1 lc rgb "#d0d0d0"

set key top right box opaque font ",14"

set xlabel "Time" font ",18"
set ylabel "Energy Variation" font ",18"

# Uncomment if desired
# set format y "%.3e"

#---------------------------------------------------------
# Line styles
#---------------------------------------------------------
set style line 1 lc rgb "#1f77b4" lw 2.5 dt 1 pt 7  ps 0.9
set style line 2 lc rgb "#ff7f0e" lw 2.5 dt 2 pt 5  ps 0.9
set style line 3 lc rgb "#2ca02c" lw 2.5 dt 3 pt 9  ps 0.9
set style line 4 lc rgb "#d62728" lw 2.5 dt 4 pt 13 ps 0.9
set style line 5 lc rgb "black"   lw 4.0 dt 1

#---------------------------------------------------------
# Initial totalenergyes (first data row = line 3)
#---------------------------------------------------------
m0 = real(system("awk 'NR==3{print $3}' results/totalenergy_0.txt"))
m1 = real(system("awk 'NR==3{print $3}' results/totalenergy_1.txt"))
m2 = real(system("awk 'NR==3{print $3}' results/totalenergy_2.txt"))
m3 = real(system("awk 'NR==3{print $3}' results/totalenergy_3.txt"))

mtot0 = m0 + m1 + m2 + m3

#---------------------------------------------------------
# Plot
#---------------------------------------------------------
plot \
    "results/totalenergy_0.txt" \
        using 2:($3-m0) every ::2 \
        with linespoints ls 1 pointinterval 20 \
        title "Zone 0", \
    "results/totalenergy_1.txt" \
        using 2:($3-m1) every ::2 \
        with linespoints ls 2 pointinterval 20 \
        title "Zone 1", \
    "results/totalenergy_2.txt" \
        using 2:($3-m2) every ::2 \
        with linespoints ls 3 pointinterval 20 \
        title "Zone 2", \
    "results/totalenergy_3.txt" \
        using 2:($3-m3) every ::2 \
        with linespoints ls 4 pointinterval 20 \
        title "Zone 3", \
    "<paste results/totalenergy_0.txt results/totalenergy_1.txt results/totalenergy_2.txt results/totalenergy_3.txt | tail -n +3" \
        using 2:(($3+$8+$13+$18)-mtot0) \
        with lines ls 5 \
        title "Total Energy"

unset output
