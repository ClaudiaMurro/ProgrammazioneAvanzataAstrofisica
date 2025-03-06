/*
 * Simulazione 1D di un tubo di shock con il metodo Smoothed Particle Hydrodynamics (SPH)
 * --------------------------------------------------------------------------------------
 * Questo codice implementa una simulazione 1D di un tubo di shock utilizzando il metodo SPH
 * con kernel cubic spline. Il dominio è periodico e le condizioni iniziali sono 
 * una regione ad alta pressione e densità a sinistra 
 * e una regione a bassa pressione e densità a destra.
 *
 * Le fasi principali della simulazione sono:
 *  1. Inizializzazione delle condizioni iniziali (densità, pressione, velocità, energia interna).
 *  2. Identificazione dei vicini per ciascuna particella utilizzando un criterio di distanza (con un qsort).
 *  3. Calcolo della densità con il metodo SPH.
 *  4. Calcolo della pressione utilizzando un'equazione di stato.
 *  5. Calcolo dell'accelerazione e della variazione di energia interna.
 *  6. Evoluzione temporale tramite il metodo di integrazione Kick-Drift-Kick (KDK).
 *  7. Salvataggio dei dati delle particelle e dei vicini a intervalli regolari.
 * Il passo temporale è determinato dal criterio di Courant.
 * La simulazione termina quando viene raggiunto il tempo finale specificato.
 *
 * Compilo: gcc sph.c -lm
 * Eseguo:  ./a.out
 * Grafico: gnuplot grafici.gp animazione.gp 
 *
 * Autore: Claudia Murro
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 3600                        // Numero totale di particelle
#define MASS 0.00015625               // Massa di ciascuna particella [kg]
#define GAMMA (5.0/3.0)               // gas monoatomico ideale
#define LZ_DENSITY 1.0                // Densità iniziale nella zona sinistra (LZ) [kg/m³]
#define RZ_DENSITY 0.125              // Densità iniziale nella zona destra (RZ) [kg/m³]
#define LZ_PRESSURE 1.0               // Pressione iniziale nella zona sinistra (LZ) (unità normalizzate)
#define RZ_PRESSURE 0.1               // Pressione iniziale nella zona destra (RZ) (unità normalizzate)
#define LZ_VELOCITY 0.0               // Velocità iniziale nella zona sinistra (LZ) [m/s]
#define RZ_VELOCITY 0.0               // Velocità iniziale nella zona destra (RZ) [m/s]
#define TUBE_LENGTH 1.0               // Lunghezza totale del tubo [m]
#define EPSILON 0.01                  // Costante per evitare problemi numerici
#define ALPHA 1.0                     // Coefficiente della viscosità artificiale
#define BETA 2.0                      // Coefficiente della viscosità artificiale
#define PRESSURE_CONVERSION 100000.0  // Conversione della pressione in Pascal
#define ENERGY_CONVERSION 100000.0    // Conversione dell'energia in Joule per chilogrammo
#define SNAPSHOT_TIME 0.01            // Tempo a cui vengono salvati gli snapshot [s]
#define END_TIME 0.08                 // Tempo finale della simulazione [s]
#define K 0.1                         // Costante per la dimensione del passo temporale
#define MAX_NEIGHBORS 64              // Numero di vicini per particella

// Struttura che rappresenta una particella SPH
typedef struct {
    double x;        // Posizione lungo il tubo 
    double rho;      // Densità 
    double velocity; // Velocità 
    double pressure; // Pressione 
    double accel;    // Accelerazione 
    double u;        // Energia interna specifica 
    double du_dt;    // Derivata temporale dell'energia interna
    double h;        // Lunghezza di smoothing 
} Particle;

// Struttura per memorizzare la lista dei vicini di una particella
typedef struct {
    int count;                      		     // Numero di vicini 
    int neighbors[MAX_NEIGHBORS];    // Indici delle particelle vicine
} NeighborList;

// Dichiarazione degli array globali per le particelle e le loro liste di vicini
Particle particles[N];               // Array delle particelle
NeighborList neighbor_lists[N];       // Array delle liste di vicini per ogni particella

//__________________________________________________________________________________KERNEL     
// Funzione segno: restituisce +1 se x > 0, -1 se x < 0, 0 se x = 0
double sgn(double x) {
    return (x > 0) - (x < 0);  // Determina il segno di x
}

// Kernel Cubic Spline per SPH
// L'argomento x rappresenta la differenza tra le posizioni di due particelle (x_i - x_j)
double cubic_spline(double x, double h) {
    double r = fabs(x);  // Consideriamo solo la distanza assoluta tra le particelle
    double q = r / h;    // Distanza normalizzata rispetto a h
    double norm = 2.0 / (3.0 * h);  // Fattore di normalizzazione

    if (q <= 1.0) {
        return norm * (1.0 - 1.5 * q * q + 0.75 * q * q * q);  // Primo intervallo del kernel
    } else if (q <= 2.0) {
        double term = 2.0 - q;
        return norm * 0.25 * term * term * term;  // Secondo intervallo del kernel
    }
    return 0.0;  // Il kernel è nullo oltre 2h
}

// Derivata del Kernel Cubic Spline
// L'argomento x rappresenta la differenza tra le posizioni di due particelle (x_i - x_j)
double grad_cubic_spline(double x, double h) {
    double r = fabs(x);  // Consideriamo solo la distanza assoluta tra le particelle
    double q = r / h;    // Distanza normalizzata rispetto a h
    double norm = 2.0 / (3.0 * h * h);  // Fattore di normalizzazione della derivata
    double grad;

    if (q <= 1.0) {
        grad = norm * (-3.0 * q + 2.25 * q * q);  // Primo intervallo della derivata
    } else if (q <= 2.0) {
        double term = 2.0 - q;
        grad = -norm * 0.75 * term * term;  // Secondo intervallo della derivata
    } else {
        grad = 0.0;  // La derivata è nulla oltre 2h
    }

    grad *= sgn(x);  // Riporta il segno corretto in base alla direzione della distanza

    return grad;
}

//__________________________________________________________________________________scrivo vicini e particelle in file, implemento c.i.

// Funzione per scrivere su file la lista dei vicini di ogni particella
void write_neighbors(int step_counter) {
    char filename[50];
    sprintf(filename, "neighbors_step_%d.dat", step_counter); // Genera nome file con il numero del passo corrente
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Errore nell'aprire il file\n");
        return;
    }

    // Scrive per ogni particella la lista dei suoi vicini
    for (int i = 0; i < N; i++) {
        fprintf(file, "Particella %d: ", i);
        for (int j = 0; j < neighbor_lists[i].count; j++) {
            fprintf(file, "%d ", neighbor_lists[i].neighbors[j]); // Stampa gli indici dei vicini
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

// Funzione per impostare le condizioni iniziali dello shock tube
void setup_shocktube() {
    // Determina il numero di particelle nelle due zone (LZ = zona sinistra, RZ = zona destra)
    int N_LZ = (int)(N * (LZ_DENSITY / (LZ_DENSITY + RZ_DENSITY)));  
    int N_RZ = N - N_LZ;  // Il resto delle particelle va nella RZ
    printf("Numero particelle LZ e RZ: %d, %d\n", N_LZ, N_RZ);

    // Inizializza le particelle nella zona di sinistra (LZ)
    for (int i = 0; i < N_LZ; i++) {
        particles[i] = (Particle){
            .x = i * (TUBE_LENGTH / (2 * N_LZ)),  // Posizione iniziale
            .rho = LZ_DENSITY,                   // Densità iniziale
            .velocity = LZ_VELOCITY,             // Velocità iniziale
            .pressure = LZ_PRESSURE,             // Pressione iniziale
            .u = LZ_PRESSURE / (LZ_DENSITY * (GAMMA - 1))  // Energia interna
        };
    }

    // Inizializza le particelle nella zona di destra (RZ)
    for (int i = N_LZ; i < N; i++) {
        particles[i] = (Particle){
            .x = (i - N_LZ) * (TUBE_LENGTH / (2 * N_RZ)) + (TUBE_LENGTH / 2),  // Posizione iniziale
            .rho = RZ_DENSITY,                   // Densità iniziale
            .velocity = RZ_VELOCITY,             // Velocità iniziale
            .pressure = RZ_PRESSURE,             // Pressione iniziale
            .u = RZ_PRESSURE / (RZ_DENSITY * (GAMMA - 1))  // Energia interna
        };
    }
}

// Funzione per scrivere i dati delle particelle su file
void write_particles(int step_counter) {
    char filename[50];
    sprintf(filename, "particles_step_%d.dat", step_counter); // Genera nome file con il numero del passo corrente
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        fprintf(stderr, "Errore: impossibile aprire il file %s per la scrittura.\n", filename);
        return;
    }

    // Scrive la riga di intestazione con le unità di misura
    fprintf(file, "# x [m]   rho [kg/m^3]   velocity [m/s]   pressure [Pa]   accel [m/s^2]   energy [J/kg]   du_dt [J/(kg·s)]   h [m]\n");

    // Scrive i dati di ogni particella
    for (int i = 0; i < N; i++) {
        double converted_pressure = particles[i].pressure * PRESSURE_CONVERSION; // Converte la pressione in Pascal
        double converted_energy = particles[i].u * ENERGY_CONVERSION; // Converte l'energia in Joule/kg
        fprintf(file, "%f %f %f %f %f %f %f %f\n", 
                particles[i].x, 
                particles[i].rho, 
                particles[i].velocity, 
                converted_pressure, 
                particles[i].accel, 
                converted_energy, 
                particles[i].du_dt,
                particles[i].h);
    }

    fclose(file);
}

//__________________________________________________________________________________calcolo densità, pressione, cs, viscosità, du/dt, dv/dt

// Applica le condizioni periodiche al valore della distanza dx tra due particelle
double apply_periodic_boundary(double dx) {
    if (dx > 0.5 * TUBE_LENGTH) {  
        dx -= TUBE_LENGTH;  // Se dx è maggiore di metà tubo, la particella è oltre il bordo destro allora riporta a sinistra
    } else if (dx < -0.5 * TUBE_LENGTH) {  
        dx += TUBE_LENGTH;  // Se dx è minore di -metà tubo, la particella è oltre il bordo sinistro allora riporta a destra
    }
    return dx;
}

// Funzione per calcolare la densità di ogni particella usando il kernel SPH
void compute_densities() {
    for (int i = 0; i < N; i++) {
        particles[i].rho = 0.0;  // Inizializza la densità della particella i
        for (int j = 0; j < neighbor_lists[i].count; j++) {
            int b = neighbor_lists[i].neighbors[j];  // Ottiene l'indice della particella vicina
            double dx = particles[i].x - particles[b].x;  // Calcola la distanza lungo x
            dx = apply_periodic_boundary(dx);  // Corregge dx per le condizioni periodiche
            
            // Determina lo smoothing length da usare
            double H = particles[i].h;  // Alternativa: usare la media degli smoothing lengths: 0.5 * (particles[i].h + particles[b].h)
            
            // Aggiorna la densità sommando il contributo della particella b pesato dal kernel
            particles[i].rho += MASS * cubic_spline(dx, H);
        }
    }
}

// Funzione per calcolare la pressione di ogni particella usando l'equazione di stato
void compute_pressures() {
    for (int i = 0; i < N; i++) {
        particles[i].pressure = (GAMMA - 1) * particles[i].rho * particles[i].u;  
        // Usa l'equazione di stato per gas perfetto
    }
}

// Calcola la velocità del suono per una particella i
double sound_speed(int i) {
    // Verifica che la pressione sia valida
    if (particles[i].pressure <= 0) {
        printf("Attenzione: pressione non valida per la particella %d (pressione = %f)\n", i, particles[i].pressure);
        return 0;  // Restituisce 0 in caso di errore di pressione
    }
    
    // Verifica che la densità sia valida
    if (particles[i].rho <= 0) {
        printf("Attenzione: densità non valida per la particella %d (densità = %f)\n", i, particles[i].rho);
        return 0;  // Restituisce 0 in caso di errore di densità
    }

    double cs = sqrt(GAMMA * particles[i].pressure / particles[i].rho);
    return cs;
}

// Calcola la viscosità artificiale tra due particelle a e b
double artificial_viscosity(int a, int b) {
    // Calcola la differenza di velocità tra le due particelle
    double v_ab = particles[a].velocity - particles[b].velocity;
    
    // Calcola la distanza tra le particelle
    double dx = particles[a].x - particles[b].x;
    dx = apply_periodic_boundary(dx);  // Applica le condizioni periodiche

    // Calcola la velocità del suono media tra le due particelle
    double cs_ab = (sound_speed(a) + sound_speed(b)) / 2.0;
    
    // Sceglie lo smoothing length per la particella a (alternativamente usare la media tra a e b)
    double H = particles[a].h;  

    if (v_ab * dx < 0) {
        double mu_ab = H * v_ab * dx / (dx * dx + EPSILON * H * H);  
        double rho_ab = (particles[a].rho + particles[b].rho) / 2.0;  // Media della densità

        // Restituisce la viscosità artificiale
        return (-ALPHA * mu_ab * cs_ab + BETA * mu_ab * mu_ab) / rho_ab;
    }

    return 0.0;  
}

// Derivata di u (energia interna) per ciascuna particella
void compute_energy_derivative() {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;  // Inizializza la somma per la derivata dell'energia
        for (int j = 0; j < neighbor_lists[i].count; j++) {
            int b = neighbor_lists[i].neighbors[j];  // Ottieni l'indice della particella vicina

            // Calcola la distanza tra le particelle con condizioni periodiche
            double dx = particles[i].x - particles[b].x;
            dx = apply_periodic_boundary(dx);  

            double H = particles[i].h;  // Smoothing length della particella i (o media tra i due smoothing lengths)
            double gradW = grad_cubic_spline(dx, H);  // Calcola il gradiente del kernel

            // Calcola la viscosità artificiale tra le particelle
            double pi_ab = artificial_viscosity(i, b);

            // Calcola la differenza di velocità tra le particelle
            double v_ab = particles[i].velocity - particles[b].velocity;

            sum += MASS * (particles[i].pressure / (particles[i].rho * particles[i].rho) + pi_ab / 2.0) * v_ab * gradW;
        }
        particles[i].du_dt = sum;
    }
}

// Calcolo dell'accelerazione per ciascuna particella
void compute_accel() {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;  // Inizializza la somma per l'accelerazione
        for (int j = 0; j < neighbor_lists[i].count; j++) {
            int b = neighbor_lists[i].neighbors[j];  // Ottieni l'indice della particella vicina
            
            // Calcola la distanza tra le particelle con condizioni periodiche
            double dx = particles[i].x - particles[b].x;
            dx = apply_periodic_boundary(dx);  

            double H = particles[i].h;  // Smoothing length della particella i (o media tra i due smoothing lengths)
            double gradW = grad_cubic_spline(dx, H);  // Calcola il gradiente del kernel

            // Calcola la viscosità artificiale tra le particelle
            double pi_ab = artificial_viscosity(i, b);

            double term = (particles[i].pressure / (particles[i].rho * particles[i].rho) + 
                           particles[b].pressure / (particles[b].rho * particles[b].rho) + pi_ab);
            
            sum -= MASS * term * gradW;
        }
        particles[i].accel = sum;
    }
}

//__________________________________________________________________________________ricerca vicini
/* Ricerca dei vicini per ciascuna particella:
 * 1. Per ogni particella i, si calcolano le distanze da tutte le altre particelle j.
 * 2. Le distanze vengono corrette applicando condizioni periodiche al dominio.
 * 3. Si memorizzano gli indici delle particelle e le rispettive distanze in un array temporaneo.
 * 4. L'array viene ordinato in base alla distanza crescente utilizzando qsort.
 * 5. Si selezionano esattamente 64 vicini più prossimi per ogni particella.
 * 6. Il valore di h di ogni particella viene aggiornato come metà della distanza del 64° vicino.
 */
 
// Struttura per ordinare le particelle in base alla distanza
typedef struct {
    int index;         // Indice della particella
    double distance;   // Distanza dalla particella di riferimento
} NeighborDistance;

// Funzione di confronto per qsort (ordina in base alla distanza crescente)
int compare_distances(const void *a, const void *b) {
    double d1 = ((NeighborDistance *)a)->distance;  // Distanza della prima particella
    double d2 = ((NeighborDistance *)b)->distance;  // Distanza della seconda particella
    return (d1 > d2) - (d1 < d2);  // Restituisce -1, 0 o 1 per il confronto
}

// Funzione per trovare i vicini di ciascuna particella con il parametro h
void find_neighbors_with_h() {
    NeighborDistance distances[N]; // Array temporaneo per memorizzare le distanze e gli indici

    for (int i = 0; i < N; i++) {
        double x_i = particles[i].x;  
        int num_neighbors = 0;

        // Calcola le distanze con tutte le altre particelle, considerando le condizioni periodiche
        for (int j = 0; j < N; j++) {
            if (i == j) continue;  // Salta se stesso

            double dx = particles[j].x - x_i;  // Distanza tra le particelle i e j
            dx = apply_periodic_boundary(dx);  // Applica la condizione di periodicità

            // Memorizza l'indice e la distanza nella struttura temporanea
            distances[num_neighbors].index = j;
            distances[num_neighbors].distance = fabs(dx);
            num_neighbors++;  // Incrementa il numero di vicini trovati
        }

        // Ordina le distanze in ordine crescente usando qsort
        qsort(distances, num_neighbors, sizeof(NeighborDistance), compare_distances);

        // Se ci sono meno di 64 vicini, stampa un errore e termina il programma
        if (num_neighbors < MAX_NEIGHBORS) {
            printf("Errore: Particella %d ha solo %d vicini!\n", i, num_neighbors);
            exit(1);  // Uscita per errore
        }

        // Seleziona esattamente 64 vicini (il massimo definito)
        for (int j = 0; j < MAX_NEIGHBORS; j++) {
            neighbor_lists[i].neighbors[j] = distances[j].index;  // Memorizza l'indice del vicino
        }
        neighbor_lists[i].count = MAX_NEIGHBORS;  // Imposta il numero di vicini

        // Imposta h in modo che la distanza massima tra vicini sia 2h
        particles[i].h = 0.5 * distances[MAX_NEIGHBORS - 1].distance;
    }
}

//__________________________________________________________________________________KDK
// Variabile per il tempo corrente
double current_time = 0.0;

// Funzione per calcolare il passo temporale
double compute_delta_t() {
    // Trova il massimo valore di h tra tutte le particelle
    double H = particles[0].h;  
    for (int i = 1; i < N; i++) {
        if (particles[i].h > H) {
            H = particles[i].h;
        }
    }

    // Inizializza min_delta_t con il valore della prima particella
    double min_delta_t = K * H / sound_speed(0);

    // Trova il passo temporale minimo per tutte le particelle
    for (int i = 1; i < N; i++) {
        double cs = sound_speed(i);  // Velocità del suono per la particella i
        double delta_t = K * H / cs;  // Calcola il passo temporale per la particella i
        
        // Aggiorna min_delta_t se il passo temporale corrente è minore
        if (delta_t < min_delta_t) {
            min_delta_t = delta_t;
        }
    }
    
    // Restituisce il passo temporale minimo
    return min_delta_t;
}

// Funzione per eseguire un passo temporale secondo il metodo KDK 
void kdk_step() {
    // Calcolare il passo temporale per questa iterazione
    double delta_t = compute_delta_t();

    // K1: Calcola velocità ed energia a metà passo usando le derivate precedenti
    for (int i = 0; i < N; i++) {
        particles[i].velocity += particles[i].accel * delta_t / 2.0;  // Aggiorna la velocità
        particles[i].u += particles[i].du_dt * delta_t / 2.0;         // Aggiorna l'energia interna
    }

    // D1: Aggiorna le posizioni a t + delta_t usando le velocità calcolate
    for (int i = 0; i < N; i++) {
        particles[i].x += particles[i].velocity * delta_t;  

        // Condizioni al contorno periodiche per mantenere le particelle nel tubo
        if (particles[i].x < 0.0) {
            particles[i].x += TUBE_LENGTH; // Riavvolgi a destra
        } else if (particles[i].x >= TUBE_LENGTH) {
            particles[i].x -= TUBE_LENGTH; // Riavvolgi a sinistra
        }
    }

    // Ricalcola i vicini, densità, pressioni, accelerazioni e derivate di energia
    find_neighbors_with_h();    // Trova i vicini 
    compute_densities();       // Calcola la densità per ogni particella
    compute_pressures();       // Calcola la pressione per ogni particella
    compute_accel();           // Calcola l'accelerazione per ogni particella
    compute_energy_derivative(); // Calcola la derivata dell'energia interna

    // K2: Calcola velocità ed energia a t + delta_t/2 con le nuove derivate
    for (int i = 0; i < N; i++) {
        particles[i].velocity += particles[i].accel * delta_t / 2.0;  // Aggiorna la velocità
        particles[i].u += particles[i].du_dt * delta_t / 2.0;         // Aggiorna l'energia interna
    }

    // Aggiorna il tempo corrente con il passo temporale
    current_time += delta_t;
}

//__________________________________________________________________________________MAIN
int main() {
    // Inizializzazione del tempo per gli snapshot
    double snapshot_time = SNAPSHOT_TIME;

    // Contatore dei passi temporali
    int step_counter = 0;

    // Imposta le condizioni iniziali di shock tube
    setup_shocktube();
    
    // Trova i vicini per ciascuna particella 
    find_neighbors_with_h();
    
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
