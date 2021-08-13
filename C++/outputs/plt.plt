set terminal png
set style line lc "orange"
set output './images/Gamma-Lyap-Exp.png'
set xlabel "Gamma"
set ylabel "Lyapunov Exponent"
plot './txt/Gamma-Lyapunov-exponents.dat' using 1:2 w lines  title "X", './txt/Gamma-Lyapunov-exponents.dat' using 1:3 w lines  title "Y",'./txt/Gamma-Lyapunov-exponents.dat' using 1:4 w lines  title "Px", './txt/Gamma-Lyapunov-exponents.dat' using 1:5 w lines title "Py"
set output './images/Gamma-Lyap-Num.png'
set xlabel "Gamma"
set ylabel "Lyapunov Numbers"
plot './txt/Gamma-Lyapunov-numbers.dat' using 1:2 w lines  title "X", './txt/Gamma-Lyapunov-numbers.dat' using 1:3 w lines  title "Y",'./txt/Gamma-Lyapunov-numbers.dat' using 1:4 w lines  title "Px", './txt/Gamma-Lyapunov-numbers.dat' using 1:5 w lines title "Py"
set output './images/Gamma-Lyap-Sum.png'
set xlabel "Gamma"
set ylabel "Lyapunov exponents sum"
plot './txt/Gamma-Lyapunov-sum.dat' w lines title "Sum of Lyapunov exponents"