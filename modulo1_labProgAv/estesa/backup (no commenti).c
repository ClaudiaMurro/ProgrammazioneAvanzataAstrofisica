#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_ATOMI 1000  // Numero totale di atomi 
#define SEMIVITA 4.5  // Tempo di dimezzamento in Gyr
#define M_SIMULAZIONI 5000  // Numero di simulazioni da eseguire
#define SEED 12345  // Seed per la generazione di numeri casuali
#define N_TEORICO 500  // Numero di punti teorici da generare

// Funzione per generare un numero casuale tra 0 e 1
double random_double() {
    return (double)rand() / RAND_MAX;
}

// Funzione per generare un numero casuale con distribuzione esponenziale di parametro lambda
double random_exponential(double lambda) {
    return -log(1.0 - random_double()) / lambda;
}

// Struttura per immagazzinare la media e la deviazione standard dei risultati delle simulazioni
typedef struct {
    double media;
    double deviazione_standard;
} Risultati;

// Funzione che simula il decadimento usando il primo metodo (probabilità di decadere)
Risultati simulazione_decadimento_classico(double t) {
    double lambda = log(2) / SEMIVITA;
    double somma = 0.0, somma_quadrati = 0.0;

    for (int m = 0; m < M_SIMULAZIONI; m++) {
        int n_decaduti = 0;
        for (int i = 0; i < N_ATOMI; i++) {
            if (random_double() < 1 - exp(-lambda * t)) {
                n_decaduti++;
            }
        }
        somma += n_decaduti;
        somma_quadrati += n_decaduti * n_decaduti;
    }

    double media = somma / M_SIMULAZIONI;
    double varianza = (somma_quadrati / M_SIMULAZIONI) - (media * media);
    double deviazione_standard = sqrt(varianza);

    return (Risultati){media, deviazione_standard};
}

// Funzione che simula il decadimento usando il secondo metodo (tempi di vita esponenziali)
Risultati simulazione_decadimento_esponenziale(double t) {
    double lambda = log(2) / SEMIVITA;
    double somma = 0.0, somma_quadrati = 0.0;

    for (int m = 0; m < M_SIMULAZIONI; m++) {
        int n_decaduti = 0;
        for (int i = 0; i < N_ATOMI; i++) {
            double tempo_decadimento = random_exponential(lambda);
            if (tempo_decadimento <= t) {
                n_decaduti++;
            }
        }
        somma += n_decaduti;
        somma_quadrati += n_decaduti * n_decaduti;
    }

    double media = somma / M_SIMULAZIONI;
    double varianza = (somma_quadrati / M_SIMULAZIONI) - (media * media);
    double deviazione_standard = sqrt(varianza);

    return (Risultati){media, deviazione_standard};
}

// Funzione per scrivere i risultati su un file .dat
void scrivi_file(char* nome_file, double tempi[], int num_tempi, Risultati risultati[], const char* tipo_simulazione) {
    FILE* file = fopen(nome_file, "w");
    if (file == NULL) {
        printf("Errore nell'aprire il file %s\n", nome_file);
        return;
    }

    fprintf(file, "# Tempo, Media, Deviazione Standard\n");
    for (int i = 0; i < num_tempi; i++) {
        fprintf(file, "%.1f %.2f %.2f\n", tempi[i], risultati[i].media, risultati[i].deviazione_standard);
    }

    fclose(file);
    printf("I risultati della %s sono stati salvati in %s\n", tipo_simulazione, nome_file);
}

// Funzione per generare i risultati teorici per N_TEORICO punti
void scrivi_file_teorico(char* nome_file, int num_punti) {
    FILE* file_teorico = fopen(nome_file, "w");
    if (file_teorico == NULL) {
        printf("Errore nell'aprire il file %s\n", nome_file);
        return;
    }

    fprintf(file_teorico, "# Tempo, Teorico\n");

    // Scriviamo la curva teorica nel file teorico
    for (int i = 0; i < num_punti; i++) {
        double t = (double)i / (num_punti - 1) * 10.0; // Calcoliamo un tempo nell'intervallo [0, 10] Gyr
        double teorico = N_ATOMI * (1 - exp(-log(2) / SEMIVITA * t));
        fprintf(file_teorico, "%.2f %.2f\n", t, teorico);
    }

    fclose(file_teorico);
    printf("I risultati teorici sono stati salvati in %s\n", nome_file);
}

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
    printf("I file con i dati sono stati generati.\n");
    return 0;
}
