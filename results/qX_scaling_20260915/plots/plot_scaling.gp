if (!exists("DATA"))   DATA   = "scaling.dat"
if (!exists("OUTPDF")) OUTPDF = "qX_Nv100_scaling.pdf"
if (!exists("OUTPNG")) OUTPNG = "qX_Nv100_scaling.png"

set title "qX p2n100: Nv=100 process-grid scaling"
set xlabel "compute nodes"
set ylabel "wall time [s]"

set xrange [3:17]
set xtics (4,8,16)
set grid
set key top right

set terminal pdfcairo enhanced size 7.5in,4.8in
set output OUTPDF

plot DATA using 1:2 with linespoints pt 7 ps 0.8 lw 1.5 \
     title "measured", \
     DATA using 1:3 with linespoints pt 5 ps 0.7 lw 1.2 dashtype 2 \
     title "ideal scaling from 4 nodes"

set terminal pngcairo enhanced size 1500,960
set output OUTPNG

replot

unset output

print sprintf("PDF: %s",OUTPDF)
print sprintf("PNG: %s",OUTPNG)
