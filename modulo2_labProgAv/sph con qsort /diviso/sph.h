#ifndef SPH
#define SPH

// Definizioni globali
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

// Strutture
typedef struct {
    double x;          // Posizione della particella
    double velocity;   // Velocità della particella
    double rho;        // Densità della particella
    double pressure;   // Pressione della particella
    double u;          // Energia interna della particella
    double accel;      // Accelerazione della particella
    double du_dt;      // Derivata dell'energia interna
    double h;          // Smoothing length della particella
} Particle;

typedef struct {
    int neighbors[MAX_NEIGHBORS]; // Indici dei vicini
    int count;                    // Numero di vicini
} NeighborList;

// Dichiarazione variabili globali
extern Particle particles[N];        // Array di particelle
extern NeighborList neighbor_lists[N];  // Lista di vicini per ogni particella
extern double current_time;          // Tempo corrente della simulazione

// Funzioni dichiarate
double sgn(double x);
void setup_shocktube();
double cubic_spline(double x, double h);
double grad_cubic_spline(double x, double h);

void find_neighbors_with_h();
void write_neighbors(int step_counter);
void write_particles(int step_counter);
double apply_periodic_boundary(double dx);
void compute_densities();
void compute_pressures();
void compute_accel();
void compute_energy_derivative();
double sound_speed(int i);
double artificial_viscosity(int a, int b);

int compare_distances(const void *a, const void *b);
double compute_delta_t();
void kdk_step();

#endif // SPH
