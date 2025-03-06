# Impostazioni generali
set terminal pngcairo enhanced size 1000,800
set output 'risultato.png'

# Impostazioni dei limiti dell'asse x
set xrange [0.3:0.8]

# Layout dei grafici
set multiplot layout 2,2 title "Risultato simulazione" font ",14"

# Grafico 1: Densità vs Posizione
set title "Densità vs Posizione"
set xlabel "Posizione (x)"
set ylabel "Densità (kg/m³)"
plot "particles_step_8.dat" using 1:2 with points pt 7 ps 1 title 'Densità'

# Grafico 2: Velocità vs Posizione
set title "Velocità vs Posizione"
set xlabel "Posizione (x)"
set ylabel "Velocità (m/s)"
plot "particles_step_8.dat" using 1:3 with points pt 7 ps 1 title 'Velocità'

# Grafico 3: Pressione vs Posizione
set title "Pressione vs Posizione"
set xlabel "Posizione (x)"
set ylabel "Pressione (Pa)"
plot "particles_step_8.dat" using 1:4 with points pt 7 ps 1 title 'Pressione'

# Grafico 4: Energia interna vs Posizione
set title "Energia interna vs Posizione"
set xlabel "Posizione (x)"
set ylabel "Energia interna (J/kg)"
plot "particles_step_8.dat" using 1:6 with points pt 7 ps 1 title 'Energia interna'

unset multiplot
