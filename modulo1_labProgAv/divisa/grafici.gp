set terminal pngcairo enhanced font 'Arial,12'

# Primo grafico: Metodo classico vs Modello teorico
set output 'Metodo1_vs_SoluzioneAnalitica.png'
set title 'Decadimento Radioattivo - Metodo 1 vs Soluzione Analitica'
set xlabel 'Tempo (Gyr)'
set ylabel 'Numero di atomi decaduti'
set grid
set key top left  # Sposta la legenda in alto a sinistra
plot 'risultati_classico.dat' using 1:2:3 with yerrorbars title 'Simulazione Metodo 1', \
     'risultati_teorici.dat' using 1:2 with lines linewidth 2 title 'Modello Teorico', \
     '-' with points pt 7 ps 1 lc rgb "red" title "Dimezzamento"
4.5 500
 e

# Secondo grafico: Metodo esponenziale vs Modello teorico
set output 'Metodo2_vs_SoluzioneAnalitica.png'
set title 'Decadimento Radioattivo - Metodo 2 vs Soluzione Analitica'
set key top left  # Sposta la legenda in alto a sinistra
plot 'risultati_esponenziale.dat' using 1:2:3 with yerrorbars title 'Simulazione Metodo 2', \
     'risultati_teorici.dat' using 1:2 with lines linewidth 2 title 'Modello Teorico', \
     '-' with points pt 7 ps 1 lc rgb "red" title "Dimezzamento"
4.5 500
 e
