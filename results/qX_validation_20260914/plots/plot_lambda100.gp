if (!exists("DATA")) DATA = "lambda100_PRODUCTION_qX_p2n100_Nv100.txt"
if (!exists("OUTPDF")) OUTPDF = "lambda100_qX_p2n100.pdf"
if (!exists("OUTPNG")) OUTPNG = "lambda100_qX_p2n100.png"

# Determine mean in diagnostic bulk region 16 <= t <= 175
stats DATA using (($1 >= 16 && $1 <= 175) ? $2 : 1/0) nooutput
bulk_mean = STATS_mean

set xlabel "temporal slice t/a"
set ylabel "lambda_{100}"
set xrange [0:191]

set title "qX p2n100: 100th covariant-Laplacian eigenvalue"

set grid
set key top center

# Mark the diagnostic bulk region
set object 1 rectangle from 16, graph 0 to 175, graph 1 behind \
    fillstyle solid 0.08 noborder

set arrow 1 from 0,bulk_mean to 191,bulk_mean \
    nohead dashtype 2

set terminal pdfcairo enhanced size 8in,4.8in
set output OUTPDF

plot DATA using 1:2 with linespoints pt 7 ps 0.35 lw 1 \
     title "lambda_{100}(t)", \
     bulk_mean with lines dashtype 2 \
     title sprintf("bulk mean = %.8f",bulk_mean)

set terminal pngcairo enhanced size 1600,960
set output OUTPNG

replot

unset output

print sprintf("bulk mean = %.17e",bulk_mean)
print sprintf("PDF: %s",OUTPDF)
print sprintf("PNG: %s",OUTPNG)
