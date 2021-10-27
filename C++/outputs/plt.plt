set terminal png
#set style line lc "orange"
set output './images/Gamma-Lyap-Exp4.png'
set xlabel "Gamma"
set ylabel "Lyapunov Exponent"
plot './txt/Gamma-Lyapunov-exponents4.dat' using 1:2 w lines  title "X", './txt/Gamma-Lyapunov-exponents4.dat' using 1:3 w lines  title "Y",'./txt/Gamma-Lyapunov-exponents4.dat' using 1:4 w lines  title "Px", './txt/Gamma-Lyapunov-exponents4.dat' using 1:5 w lines title "Py"
set output './images/Gamma-Lyap-Num4.png'
set xlabel "Gamma"
set ylabel "Lyapunov Numbers"
plot './txt/Gamma-Lyapunov-numbers4.dat' using 1:2 w lines  title "X", './txt/Gamma-Lyapunov-numbers4.dat' using 1:3 w lines  title "Y",'./txt/Gamma-Lyapunov-numbers4.dat' using 1:4 w lines  title "Px", './txt/Gamma-Lyapunov-numbers4.dat' using 1:5 w lines title "Py"
set output './images/Gamma-Lyap-Sum4.png'
set xlabel "Gamma"
set ylabel "Lyapunov exponents sum"
plot './txt/Gamma-Lyapunov-sum4.dat' w lines title "Sum of Lyapunov exponents"