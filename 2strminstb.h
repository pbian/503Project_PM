/*   2strminstb.h Two-Stream Instability   */


// We use static arrays so no dynamic allocations
#define NMAX 10000			// for a sufficiently large array
#define mad 500				// for a small array
#define ppc 8				// particles per cell
#define tt 40				// total time steps

int Lx,dx;
int count;
int Npart;
double v_th,v_da,v_db;
double delta2;
double omega;
double dt;

int Amatrix[NMAX][NMAX] = {0};

// Two populations, A and B
// [tt] means tt time layers
// [NMAX] is number of particles
// [7] means 7 variables for each particle: index (This particle is the index_th particle in the TOTAL particles); x, y, z; Vx, Vy, Vz
double Part_Matrix_A[tt][NMAX][7];
double Part_Matrix_B[tt][NMAX][7];

// charge for populations A and B
double chargenodeA[tt][NMAX];
double chargenodeB[tt][NMAX];

// normalized charges
double normchargenodeA[tt][NMAX];
double normchargenodeB[tt][NMAX];

double normdiff[tt][NMAX];

double efield[tt][NMAX];
double phi[tt][NMAX];

void setup_domain();
void nodecharge_efield();
void finite_diff_matrix();
void periodic_move_node();
void output();
void vel_temp();
void energy();

double update_phi_field(double rho_phi[]);
double update_E_field(double phi_ef[]);
double rng(double num);
