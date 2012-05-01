/*   ps.h Plasma_Sheath   */

#define AMAX 320000				//	for a sufficiently large array
#define AMIN 1000				//	for a small array
#define A 50						//	for minimal array
#define indx 7					//	particle matrix index
#define tt 1000					//	max time steps

int Lx,Ly,dx,dy,delta2;
int imax,jmax,nmax;
int LeftBC,RightBC;
int ppc,Npart;
int newcount;
int finaltime;

int pct[tt] = {0};

double v_o,v_th;
double pcn;
double pi;
double nt,dt;
double omega;

double Part_Matrix[AMAX][indx];
double PMO[AMAX][indx];
double dpm[AMAX][indx];
double Amatrix[AMIN][AMIN] = {0};
double ion_node[A][A];
double norm_inode[A][A];
double phi[A][A];
double Ex[A][A];
double Ey[A][A];

void initial_conditions();
void inject_node_move();
void outputs();

double nl_phi_solver(double ni_phi[]);
double rng(double num);
double read_timer();
void * inject_node_move_parallel(void *thread_d);
void * inject_particles_parallel(void *thread_d);