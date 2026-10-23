set encoding utf8
set termoption noenhanced
set title "* 3rd order, 2db ripple, 5ghz cutoff, z0 = 50"
set xlabel "s"
set grid
unset logscale x 
set xrange [0.000000e+00:3.000000e-09]
unset logscale y 
set yrange [-1.118704e-01:1.632575e+00]
#set xtics 1
#set x2tics 1
#set ytics 1
#set y2tics 1
set format y "%g"
set format x "%g"
plot 'np_gptmp.data' using 1:2 with lines lw 1 title "v(2)",\
'np_gptmp.data' using 3:4 with lines lw 1 title "v(2)-v(1)/2"
