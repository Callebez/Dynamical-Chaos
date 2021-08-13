#set terminal png
unset title
tit(n) = sprintf("Trajectory with gamma = %f",n*0.01)
#do for[i=0:150]{
#    a=i/100
#    set output './outputs/images/Discrete/teste'.i.'.png'
#    set autoscale
#    plot './outputs/txt/Discrete/teste'.i.'.dat' using 1:2 w lines lc 14 title tit(i)
#}
set terminal gif animate delay 10
set output './outputs/gifs/Trajectory-Gamma2.gif'
set xrange [-40:40]
set yrange [-40:40]
do for [i=0:150]{
    plot './outputs/txt/Discrete/teste'.i.'.dat' using 1:2 w lines lc 14 title tit(i)
}
