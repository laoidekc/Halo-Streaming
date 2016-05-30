
# Gnuplot script file for plotting data in file "receives_no_calculation.dat"
#set size 0.5, 0.5
set term png size 600, 400
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set xlabel "Message number"
      set ylabel "Interval (s)"
      unset key
set output "extraction_output_nocalc.png"
      set title "Extraction rate (no calculation) erratic behaviour"
plot "extraction_output_nocalc.dat" using 1:3 with points
set output "extraction_output_1Mcalc.png"
      set title "Extraction rate (1M calculation) erratic behaviour"
plot "extraction_output_1Mcalc.dat" using 1:3 with points
set output "injection_output_nocalc.png"
      set title "Injection rate (no calculation) erratic behaviour"
plot "injection_output_nocalc.dat" using 1:3 with points
set output "injection_output_1Mcalc.png"
      set title "Injection rate (1M calculation) erratic behaviour"
plot "injection_output_1Mcalc.dat" using 1:3 with points
