#ifndef DECADIMENTO_H
#define DECADIMENTO_H

// Costanti
#define N_ATOMI 1000  // Numero totale di atomi 
#define SEMIVITA 4.5  // Tempo di dimezzamento in Gyr
#define M_SIMULAZIONI 5000  // Numero di simulazioni da eseguire
#define SEED 12345  // Seed per la generazione di numeri casuali
#define N_TEORICO 500  // Numero di punti teorici da generare

// Struttura per immagazzinare la media e la deviazione standard dei risultati delle simulazioni
typedef struct {
    double media;
    double deviazione_standard;
} Risultati;

// Funzione per generare un numero casuale tra 0 e 1
double random_double();

// Funzione per generare un numero casuale con distribuzione esponenziale di parametro lambda
double random_exponential(double lambda);

// Funzione che simula il decadimento usando il primo metodo (probabilità di decadere)
Risultati simulazione_decadimento_classico(double t);

// Funzione che simula il decadimento usando il secondo metodo (tempi di vita esponenziali)
Risultati simulazione_decadimento_esponenziale(double t);

// Funzione per scrivere i risultati su un file .dat
void scrivi_file(char* nome_file, double tempi[], int num_tempi, Risultati risultati[], const char* tipo_simulazione);

// Funzione per generare i risultati teorici per N_TEORICO punti
void scrivi_file_teorico(char* nome_file, int num_punti);

#endif
