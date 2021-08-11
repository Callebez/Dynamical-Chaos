set terminal png
do for[i=0:59]{
    set output './outputs/images/Discrete/teste'.i.'.png'
    set autoscale
    plot './outputs/txt/Discrete/teste'.i.'.dat' using 1:2 w lines
}