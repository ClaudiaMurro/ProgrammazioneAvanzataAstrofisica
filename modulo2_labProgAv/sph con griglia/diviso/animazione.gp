# Impostazioni per la generazione della GIF animata
set terminal gif animate delay 100 loop 0 size 1000,800
set output "shock_tube_evolution.gif"

do for [i=0:8] {
    set multiplot layout 2,2 title sprintf("Shock Tube Evolution - Step %d", i) font ",14"

    # Grafico 1: Densità vs Posizione
    set title "Densità vs Posizione"
    set xlabel "Posizione (x)"
    set ylabel "Densità (kg/m³)"
    plot sprintf("particles_step_%d.dat", i) using 1:2 with points pt 7 ps 1 title 'Densità'

    # Grafico 2: Velocità vs Posizione
    set title "Velocità vs Posizione"
    set xlabel "Posizione (x)"
    set ylabel "Velocità (m/s)"
    plot sprintf("particles_step_%d.dat", i) using 1:3 with points pt 7 ps 1 title 'Velocità'

    # Grafico 3: Pressione vs Posizione
    set title "Pressione vs Posizione"
    set xlabel "Posizione (x)"
    set ylabel "Pressione (Pa)"
    plot sprintf("particles_step_%d.dat", i) using 1:4 with points pt 7 ps 1 title 'Pressione'

    # Grafico 4: Energia interna vs Posizione
    set title "Energia interna vs Posizione"
    set xlabel "Posizione (x)"
    set ylabel "Energia interna (J/kg)"
    plot sprintf("particles_step_%d.dat", i) using 1:6 with points pt 7 ps 1 title 'Energia interna'

    unset multiplot
}

set output  # Chiude il file GIF
