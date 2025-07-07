plot "disp_vert" using (1000*$1):(2.0*18.7257555+$2-$3) with l lw 2 title "Verticle Width"
replot "disp.101" using (1000*$1):(2.0*18.7257555+2.0*$4) with l lw 2 title "Horizontal Width"
set grid
set xrange[0:4]
set xlabel "Time (ms)"
set yrange[-10:60]
set ylabel "Width (mm)"
set term png
set output "CylinderDeformation.png"
replot

plot "pressure_probes.txt" using (1000*$2):(1.0e-6*$3) with l lw 2 title "Sensor P1"
set grid
set xrange[0:4]
set xlabel "Time (ms)"
set yrange[0:10]
set ylabel "Pressure (MPa)"
set term png
set output "SensorPressure.png"
replot
