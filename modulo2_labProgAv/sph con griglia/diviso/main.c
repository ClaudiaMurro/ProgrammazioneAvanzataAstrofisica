#include "sph.h"
#include <stdio.h>

int main() {
    // Inizializzazione del tempo per gli snapshot
    double snapshot_time = SNAPSHOT_TIME;

    // Contatore dei passi temporali
    int step_counter = 0;

    // Imposta le condizioni iniziali di shock tube
    setup_shocktube();
    
    // Calcola gli h adattivi
    compute_smoothing_length();
    
    // Trova i vicini per ciascuna particella
    find_neighbors_with_grid();
    
    // Scrivi i dati iniziali dei vicini in un file
    write_neighbors(step_counter);
    
    // Scrivi i dati iniziali delle particelle in un file
    write_particles(step_counter);

    // Calcola densità, pressioni, accelerazioni e derivate dell'energia
    compute_densities();
    compute_pressures();
    compute_accel();
    compute_energy_derivative(); 

    // Ciclo principale della simulazione, che continua fino al tempo finale
    while (current_time < END_TIME) {
        // Esegui un passo temporale KDK
        kdk_step();

        // Scrivere i dati a intervalli di tempo specificati
        if (current_time >= snapshot_time) {
            step_counter++;  // Incrementa il contatore dei passi temporali
            write_particles(step_counter);  // Scrivi i dati delle particelle con il nome dinamico
            write_neighbors(step_counter);  // Scrivi i dati dei vicini con il nome dinamico
            
            // Aggiorna il prossimo intervallo di snapshot
            snapshot_time += SNAPSHOT_TIME;

            // Stampa il passo completato e il contatore dei passi
            printf("Passo %d completato\n", step_counter);
        }
    }
    
    // Messaggio di fine simulazione
    printf("Simulazione terminata!\n");
    
    return 0;
}