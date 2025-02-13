#include <stdio.h>
#include <stdlib.h>
#include "decadimento.h"

int main() {
    srand(SEED);  // Inizializzazione del generatore di numeri casuali

    // Definiamo un array con i tempi in cui calcolare il decadimento
    double tempi[] = {1.0, 1.5, 2.0, 4.0, 4.5, 8.0};
    int num_tempi = sizeof(tempi) / sizeof(tempi[0]);

    printf("i tempi sono in Gyr\n");
    
    // Eseguiamo la simulazione con il primo metodo (probabilità di decadimento)
    Risultati risultati_classico[num_tempi];
    printf("=== Risultati Simulazione Metodo 1 (probabilità di decadere) ===\n");
    printf("Tempo \tMedia\tDev.Std\n");
    for (int i = 0; i < num_tempi; i++) {
        risultati_classico[i] = simulazione_decadimento_classico(tempi[i]);
        printf("%.1f\t%.2f\t%.2f\n", tempi[i], risultati_classico[i].media, risultati_classico[i].deviazione_standard);
    }
    scrivi_file("risultati_classico.dat", tempi, num_tempi, risultati_classico, "Simulazione Metodo 1");

    // Eseguiamo la simulazione con il secondo metodo (tempi di vita esponenziali)
    Risultati risultati_esponenziale[num_tempi];
    printf("\n=== Risultati Simulazione Metodo 2 (tempi di vita) ===\n");
    printf("Tempo \tMedia\tDev.Std\n");
    for (int i = 0; i < num_tempi; i++) {
        risultati_esponenziale[i] = simulazione_decadimento_esponenziale(tempi[i]);
        printf("%.1f\t%.2f\t%.2f\n", tempi[i], risultati_esponenziale[i].media, risultati_esponenziale[i].deviazione_standard);
    }
    scrivi_file("risultati_esponenziale.dat", tempi, num_tempi, risultati_esponenziale, "Simulazione Metodo 2");

    // Scriviamo i risultati teorici per N_TEORICO punti
    scrivi_file_teorico("risultati_teorici.dat", N_TEORICO);

    printf("\nSimulazione completata con successo!\n");
    return 0;
}
