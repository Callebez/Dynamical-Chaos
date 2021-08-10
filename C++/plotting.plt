set terminal png
set output './outputs/images/teste.png'
set autoscale
plot './outputs/txt/teste.dat' using 1:2 w lines