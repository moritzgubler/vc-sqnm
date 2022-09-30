#!/usr/bin/gnuplot

set term pdf enhanced size 8.6cm,5.733333cm font "Times,8" fontscale 0.75

set style line 1 lt 1 lc rgb "#A00000" lw 2 pt 7 ps .75
set style line 2 lt 1 lc rgb "#00A000" lw 2 pt 11 ps .75
set style line 3 lt 1 lc rgb "#5060D0" lw 2 pt 9 ps .75
set style line 4 lt 1 lc rgb "#0000A0" lw 2 pt 8 ps .75
set style line 5 lt 1 lc rgb "#D0D000" lw 2 pt 13 ps .75
set style line 6 lt 1 lc rgb "#00D0D0" lw 2 pt 12 ps .75
set style line 7 lt 1 lc rgb "#B200B2" lw 2 pt 5 ps .75

#set terminal wxt enhanced color font 'Helvetica,12' linewidth 0.6
#set terminal postscript enhanced color font 'Helvetica,12' linewidth 0.6
set output "estimates.pdf"
#set key font ",8"
#set key horizontal bottom left width 2 maxcolumns 2
set xtics rotate by 45 offset -.4,-1.5
#set title "Lattice step size and inverse of biggest EV of hessian matrix with 64 Si atoms"
set xlabel "Iteration"
set ylabel 'Energy [Ha]'
set logscale y 10
set title 'Energy estimate'
plot "estimates.txt" u 1 title "Energy" ls 1 w lp, \
     "" u 2 title "||{/Times-Italic F}|| ^2 / {/Symbol l}_{min}" ls 2 w lp