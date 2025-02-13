/*
 * Simulazione del decadimento radioattivo di N_ATOMI atomi utilizzando due metodi:
 * 1. Metodo 1: calcolo della probabilità di decadimento per ogni atomo.
 * 2. Metodo 2: generazione di tempi di decadimento esponenziali per ogni atomo.
 *
 * Funzionalità principali:
 * 1. Calcolo del decadimento in base a probabilità di decadimento o tempi di vita esponenziali.
 * 2. Calcolo della media e deviazione standard per ogni simulazione.
 * 3. Salvataggio dei risultati delle simulazioni in file .dat 
 * 4. Generazione di dati teorici (da soluzione analitica) per il confronto con le simulazioni.
 *
 * Compilo: gcc decadimento.c -lm
 * Eseguo:  ./a.out
 * Grafico: gnuplot grafici.gp
 *
 * Autore: Claudia Murro
 */

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
    double lambda = log(2) / SEMIVITA;  // Calcolo di lambda dalla semivita
    double somma = 0.0, somma_quadrati = 0.0;

    // Eseguiamo M_SIMULAZIONI simulazioni
    for (int m = 0; m < M_SIMULAZIONI; m++) {
        int n_decaduti = 0;
        // Per ogni atomo, controlliamo se è decaduto
        for (int i = 0; i < N_ATOMI; i++) {
            // Condizione per determinare se l'atomo è decaduto in un determinato tempo
            if (random_double() < 1 - exp(-lambda * t)) {
                n_decaduti++;  // Se la probabilità è inferiore alla soglia, l'atomo è decaduto
            }
        }
        somma += n_decaduti;  // Sommiamo i decaduti per questa simulazione
        somma_quadrati += n_decaduti * n_decaduti;  // Sommiamo i quadrati dei decaduti per calcolare la varianza
    }

    // Calcoliamo la media e la deviazione standard
    double media = somma / M_SIMULAZIONI;
    double varianza = (somma_quadrati / M_SIMULAZIONI) - (media * media);
    double deviazione_standard = sqrt(varianza);

    return (Risultati){media, deviazione_standard};  // Restituiamo la media e la deviazione standard
}

// Funzione che simula il decadimento usando il secondo metodo (tempi di vita esponenziali)
Risultati simulazione_decadimento_esponenziale(double t) {
    double lambda = log(2) / SEMIVITA;  // Calcolo di lambda dalla semivita
    double somma = 0.0, somma_quadrati = 0.0;

    // Eseguiamo M_SIMULAZIONI simulazioni
    for (int m = 0; m < M_SIMULAZIONI; m++) {
        int n_decaduti = 0;
        // Per ogni atomo, generiamo un tempo di decadimento e vediamo se è prima del tempo t
        for (int i = 0; i < N_ATOMI; i++) {
            // Generiamo il tempo di decadimento esponenziale per ogni atomo
            double tempo_decadimento = random_exponential(lambda);
            // Condizione per verificare se il tempo di decadimento è inferiore o uguale a t
            if (tempo_decadimento <= t) {
                n_decaduti++;  // Se il tempo di decadimento è inferiore a t, l'atomo è decaduto
            }
        }
        somma += n_decaduti;  // Sommiamo i decaduti per questa simulazione
        somma_quadrati += n_decaduti * n_decaduti;  // Sommiamo i quadrati dei decaduti per calcolare la varianza
    }

    // Calcoliamo la media e la deviazione standard
    double media = somma / M_SIMULAZIONI;
    double varianza = (somma_quadrati / M_SIMULAZIONI) - (media * media);
    double deviazione_standard = sqrt(varianza);

    return (Risultati){media, deviazione_standard};  // Restituiamo la media e la deviazione standard
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
    printf("Per graficare i risultati dei due modelli e confrontarli con la soluzione analitica eseguire: gnuplot grafici.pg\n");
    return 0;
}
