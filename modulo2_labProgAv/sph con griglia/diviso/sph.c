#include "sph.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definizioni delle variabili globali
Particle particles[N];
NeighborList neighbor_lists[N];

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
/*
 * Funzione per calcolare i vicini di ogni particella utilizzando una griglia adattiva e un h variabile.
 * 1. Compute_smoothing_length(): 
 *    - Per ogni particella, calcola il valore ottimale di "h" (lunghezza di smorzamento) in modo adattivo, cercando di mantenere un numero costante di vicini (MAX_NEIGHBORS).
 *    - Inizialmente, viene calcolato un valore iniziale di "h_start" e, attraverso un ciclo iterativo, viene regolato affinché l'intervallo [x - 2*h, x + 2*h] contenga esattamente MAX_NEIGHBORS particelle vicine. 
 *    - Se ci sono troppo pochi vicini, h_start viene aumentato; se ci sono troppi vicini, h_start viene diminuito.
 *
 * 2. find_HMAX(): 
 *    - Trova il massimo valore di "h" tra tutte le particelle, che viene poi utilizzato per determinare la dimensione delle celle della griglia.
 *
 * 3. find_neighbors_with_grid(): 
 *    - Una volta calcolato HMAX, la funzione crea una griglia di celle di dimensione 2*HMAX per gestire la ricerca dei vicini.
 *    - Le particelle vengono assegnate alle celle della griglia in base alla loro posizione.
 *    - Successivamente, per ciascuna particella, vengono cercati i vicini nelle celle di appartenenza e nelle celle adiacenti (tenendo conto delle condizioni periodiche per il calcolo delle distanze).
 *    - I vicini vengono memorizzati in una lista di vicini per ciascuna particella.
 *    - Alla fine, la memoria utilizzata per le celle della griglia viene liberata.
 */

// Funzione per calcolare la distanza periodica
double periodic_distance(double x1, double x2) {
    double dx = fabs(x1 - x2);
    if (dx > TUBE_LENGTH / 2) dx = TUBE_LENGTH - dx;
    return dx;
}

// Funzione per trovare l'h adattivo per ogni particella
void compute_smoothing_length() {
    for (int i = 0; i < N; i++) {
        double h_start = 32 * TUBE_LENGTH / N;  // Valore iniziale per h_start
        int count = 0;
        int iter = 0;
        int max_iterations = 500;  // Limite massimo di iterazioni

        // Ciclo per trovare il giusto h_start
        while (count != MAX_NEIGHBORS && iter < max_iterations) {
            iter++;
            count = 0;
            double range = 2 * h_start;  // Intervallo [x - 2*h_start, x + 2*h_start]
            
            // Calcolo dei vicini nell'intervallo
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    double dx = periodic_distance(particles[i].x, particles[j].x);
                    if (dx <= range) {
                        count++;
                    }
                }
            }

            // Adattamento di h_start
            if (count < MAX_NEIGHBORS) {
                h_start *= 1.05;  // Aumenta h_start se ci sono troppi pochi vicini
            } else if (count > MAX_NEIGHBORS) {
                h_start *= 0.95;  // Diminuisci h_start se ci sono troppi vicini
            }
        }

        // Assegna h_start come il valore finale di h per la particella
        particles[i].h = h_start;
    }
}

// Funzione per trovare il massimo valore di h
double find_HMAX() {
    double HMAX = 0.0;
    for (int i = 0; i < N; i++) {
        if (particles[i].h > HMAX) {
            HMAX = particles[i].h;
        }
    }
    return HMAX;
}

void find_neighbors_with_grid() {
    // Calcola il massimo valore di h tra tutte le particelle (HMAX)
    double HMAX = find_HMAX();
    
    // Determina il numero di celle nella griglia
    int num_cells = (int)(TUBE_LENGTH / (2 * HMAX)) + 1;
    double cell_size = 2 * HMAX;
    
    // Struttura per la lista linkata delle celle
    typedef struct CellNode {
        int index;
        struct CellNode *next;
    } CellNode;
    
    // Array di puntatori per le celle della griglia
    CellNode *grid[num_cells];
    for (int i = 0; i < num_cells; i++) {
        grid[i] = NULL;
    }
    
    // Inserisce le particelle nelle rispettive celle
    for (int i = 0; i < N; i++) {
        int cell_index = ((int)(particles[i].x / cell_size) + num_cells) % num_cells;
        
        CellNode *new_node = (CellNode *)malloc(sizeof(CellNode));
        if (!new_node) {
            fprintf(stderr, "Errore di allocazione memoria per CellNode\n");
            exit(EXIT_FAILURE);
        }

        new_node->index = i;
        new_node->next = grid[cell_index];
        grid[cell_index] = new_node;
    }
    
    // Trova i vicini per ciascuna particella
    for (int i = 0; i < N; i++) {
        int cell_index = ((int)(particles[i].x / cell_size) + num_cells) % num_cells;
        
        neighbor_lists[i].count = 0;
        
        // Controlla la cella attuale e le due celle adiacenti
        for (int offset = -1; offset <= 1; offset++) {
            int neighbor_cell = (cell_index + offset + num_cells) % num_cells;
            
            for (CellNode *node = grid[neighbor_cell]; node != NULL; node = node->next) {
                int j = node->index;
                if (i != j) { // Evita di considerare la particella stessa
                    double dx = periodic_distance(particles[i].x, particles[j].x);
                    if (dx <= 2 * particles[i].h) {
                        if (neighbor_lists[i].count < MAX_NEIGHBORS) {
                            neighbor_lists[i].neighbors[neighbor_lists[i].count++] = j;
                        }
                    }
                }
            }
        }
    }
    
    // Libera la memoria delle liste delle celle
    for (int i = 0; i < num_cells; i++) {
        CellNode *node = grid[i];
        while (node) {
            CellNode *tmp = node;
            node = node->next;
            free(tmp);
        }
    }
}

//__________________________________________________________________________________KDK
// Variabile per il tempo corrente
double current_time = 0.0;

// Funzione per calcolare il passo temporale
double compute_delta_t() {
    // Trova il massimo valore di h tra tutte le particelle
    double H = find_HMAX(); //particles[0].h;  

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
    compute_smoothing_length();
    find_neighbors_with_grid();    // Trova i vicini 
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
