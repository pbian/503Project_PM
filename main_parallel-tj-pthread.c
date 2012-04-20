// A sample Particle-in-Cell (PIC) code for Electro-Static (ES) Plasma Dynamics Simulation
// 1-D Two-Stream Instability
// For CSCI 503 Project

/*

------------Some Physical Review------------
Plasma: charged particles, e.g., ions, electrons. Ions have positive charge, while electrons have negative charge. Their motions follows the Newton's 2nd Law:

m*a = q*E
where 'm' is mass, 'a' is dv/dt (acceleration), 'q' is charge, 'E' is electric field (E field).

------------What is the code doing?------------
In this code, we are tracing the positions of electrons (negatively charged particles). Since the mass of ions are much larger than the mass of electrons, the ions are assumed to be im-mobile (fixed at their locations). So, only electrons are studied.

The problem involved in this code is called "Two-Stream Instability", and the physical description is:

Two streams of electrons are counter-flowing with each other in x-direction (1-D problem), like the following 'figure',


Population A		  Population B
-------------->		<--------------
-------------->		<--------------
-------------->		<--------------

and the coordinate system is:

y

^
|
|
|----> x

Since there is no dependence on the y-direction, this problem is 1-D. Like mentioned above, the core part is to trace the position, and the steps/maps are:

--> Initially, their positions/velocities are given (input).
--> In order to get the new position, we need to know the new velocity (x += v * delta_t);
--> In order to know the new velocity, we need to know the acceleration (v += a * delta_t);
--> In order to know the acceleration, we need to know the E field (m * a = q * E);
--> In order to know the E field, we need to know the electric potential, Greek letter 'Phi' (E = - d(Phi)/dx);
--> In order to know the electric potential, we need to know the charge density, Greek letter rho (Laplacian of Phi = rho);
--> In order to know the charge density rho, we need to know the positions of particles (so this goes back to the beginning of the cycle);

The solution is time-serial:

	Initial particles position/velocity -->

|------>Get charge density from the positions -->
|	Solve for Phi -->
|	Take the derivative of Phi, get E -->
|	Get acceleration for each particle -->
|	Get new velocity of each particle -->
|	Get new position of each particle -->
--------Next time step

The above loop is done until 'steady state' is obtained - continuing the time series, the picture stays the same pattern.


------------How to do this in the code?------------
The simulation is based on a 1-D grid:

+-----+-----+-----+-----+-----+-----+

Where '+' stands for the grid node. Like we mentioned above, the ions are fixed - so their contribution of positive charge are directly defined at the grid node (given as the input). In the code, the Phi (potential), E (E field), and rho (charge density) are defined at the grid node, and the particles information are stored in an long array. So this method is called 'Particle-In-Cell'.

So, in the code, each step is:

--Initialize particle position/velocity: This is a large array, many rows, each row stores the position and velocity.

--Get the charge density, rho, at the grid node: All the particles contribute their charge to their near nodes. For example, if one particle's position is between node i and i+1, and the charge carried by that particle will be counted onto both node i and node i+1 - with some weight, and the weight is actually the distance (if the particle is nearer to node i, then more portion of the charge is going to be at node i). This is done for all the particles, so if we parallelize this part, there IS a data-dependency problem (the rho at each node is a sum, while this sum is from all the particles nearby).

--Solve for Phi from charge density rho. This is matrix solver problem - basically, AX = B, solve for X. In this code, the size of the matrix is based on how many nodes we have (i.e., how do we discretize the 1-D simulation domain). We can also think about this part for parallelism.

--Get E from Phi - just take the derivative.

--Get acceleration for each particle, in the m*a = q*E equation, only a is unknow now - note here, E is interpolated between two nodes. For example, if the particle's current position is at the middle of two nodes, then the E field it feels is the average of the E's of the two nodes. This is done for all the particles.

--Get new velocity and new position. This is done for all the particles.


The solution loop is:

		Initialize particle position/velocity -->

	----->	Get the charge at each grid node (electrons carry charges) -->
	|
	|	Solve the Phi field from AX = B (this is a matrix solver problem) -->
	|
	|	Get the E field from Phi, E = -d(Phi)/dx -->
	|
	-------	Update particle velocity and position ( m*a = q*E, and a = dv/dt; so v += a*dt; x += v*dt ) -->






The pseudo-flowchart is in the 'main' part.

*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2strminstb.h"	// This is the header file
#include <libc.h>
#include <omp.h>
#include <pthread.h>
#define NUM_THREADS 20
pthread_mutex_t mutexcharge = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexefield;
struct thread_data {
    int threadNum;
    int numElements;
    int numElementsPerThread;
    int currentii;
};

int main () {
	//srand(time(0));
	// ++++++++++++++++++++++++++++++++
	// Inputs
	v_th = 1; v_da = 5; v_db = -5;		// initial velocity, 'th' means 'thermal (random)', 'd' means 'drifting'
        dx = 1;			// domain length and grid size
	delta2 = dx*dx;				// grid size squared
	dt = .5/v_da;				// time step size
	omega = .75;				// a factor in solving AX = B (SOR method)
	// ================================
	
	// Setup the domain
	setup_domain();

	// Form the coefficients of AX = B when solving Phi
	finite_diff_matrix();

	// Get the charge on each grid node, and solve for Phi, and get E field
	// %%%%%%%%%%%%% There is a data-dependency when getting the charge density as I mentioned above %%%%%%%%%%%%%%%%%
	nodecharge_efield();
	double seconds = read_timer();
	// Update the particle velocity/position
	// %%%%%%%%%%%%% This is the part we need to parallelize %%%%%%%%%%%%%%%%%
	periodic_move_node();
	double elapsed = read_timer() - seconds;
	// Output checkpoint data for post-processing
	output();
	vel_temp();	// velocity and temperature (temperature is related to random velocities)
	energy();	// energy
	
	FILE *file1;
	file1 = fopen("v_th.txt", "w");
	fprintf(file1,"%f\t%f\t%f\t%d\t%d",v_th,v_da,v_db,Lx,tt);
	fclose(file1);
	printf("elapsed time (parallel): %lf seconds\n", elapsed);
	return EXIT_SUCCESS;
}


double read_timer()
{
  static int initialized = 0;
  static struct timeval start;
  struct timeval end;
  if( !initialized ) {
    gettimeofday( &start, NULL );
    initialized = 1;
  }

  gettimeofday( &end, NULL );

  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}


void setup_domain() {
	
	int i,k;
	double a;
	const double pi = 4*atan(1);
	double x_pos,v_x,v_y,v_z;
	double vrand1,vrand2;
	double th1,th2;
	
	Npart = Lx*ppc; // number of particles per population
	count = 0;
	for (i=0; i<Lx; i++) {
		for (k=0; k<ppc; k++) {

			// Initialize Population A
			x_pos = i*dx + rng(a);
			vrand1 = sqrt(-log(1 - .9999999*rng(a)));
			vrand2 = sqrt(-log(1 - .9999999*rng(a)));
			th1 = 2*pi*rng(a);	th2 = 2*pi*rng(a);
			v_x = v_da + v_th*vrand1*cos(th1);
			v_y = v_th*vrand1*sin(th1);
			v_z = v_th*vrand2*cos(th2);
			
			Part_Matrix_A[0][count][0] = count + 1;
			Part_Matrix_A[0][count][1] = x_pos;	// 1-D, only x_position is needed
			Part_Matrix_A[0][count][4] = v_x;
			Part_Matrix_A[0][count][5] = v_y;
			Part_Matrix_A[0][count][6] = v_z;
			

			// Initialize Population B
			x_pos = i*dx + rng(a);
			vrand1 = sqrt(-log(1 - .9999999*rng(a)));
			vrand2 = sqrt(-log(1 - .9999999*rng(a)));
			th1 = 2*pi*rng(a);	th2 = 2*pi*rng(a);
			v_x = v_db + v_th*vrand1*cos(th1);
			v_y = v_th*vrand1*sin(th1);
			v_z = v_th*vrand2*cos(th2);
			
			Part_Matrix_B[0][count][0] = count + 1;
			Part_Matrix_B[0][count][1] = x_pos;	// 1-D, only x_position is needed
			Part_Matrix_B[0][count][4] = v_x;
			Part_Matrix_B[0][count][5] = v_y;
			Part_Matrix_B[0][count][6] = v_z;
			
			count++;}
	}
	
	printf("\n\n====System Parameters====\n\n");
	printf("Initial V_th:\t\t%8.0f\n",v_th);
	printf("Initial V_1:\t\t%8.0f\n",v_da);
	printf("Initial V_2:\t\t%8.0f\n",v_db);
	printf("Grid Resolution:\t%8d\n",dx);
	printf("System Lengths:\t%8d\n",Lx);
	printf("Particles/Cell:\t%8d\n",2*ppc);
	printf("# of Particles:\t%8d\n",2*Npart);
	printf("Time Increment:\t%8.2f\n",dt);
	printf("Final Time Step:\t%8d\n",tt);

}

void finite_diff_matrix() {
	
	//	Construct A Matrix
	int i,j; double delta2;
	delta2 = dx*dx;
	for (i=0; i<Lx; i++) {
		if (i == 0) {Amatrix[i][i] = -2; Amatrix[i][i+1] = 1; Amatrix[i][Lx-1] = 1;}
		else if (i == Lx-1) {Amatrix[i][i] = -2; Amatrix[i][i-1] = 1; Amatrix[i][0] = 1;}
		else {Amatrix[i][i] = -2; Amatrix[i][i+1] = 1; Amatrix[i][i-1] = 1;}
	}

	
	FILE *file1;
	file1 = fopen("Amatrix.txt", "w");
	for (i=0; i<Lx; i++) {
		for (j=0; j<Lx; j++) {fprintf(file1, "%2d ",Amatrix[i][j]);}
		fprintf(file1, "\n");}
	fclose(file1);
	
}

// Get the charge on each grid node, solve for E field
void nodecharge_efield() {
	
	int i,j,k,m,n,ic,ef;
	double Area1 = 0, Area2 = 0, Area3 = 0, Area4 = 0;
	double rho[Lx+1],pef[Lx+1];
	for (k=0; k<Npart; k++) {
		i = floor(Part_Matrix_A[0][k][1]); m = ceil(Part_Matrix_A[0][k][1]);
		j = floor(Part_Matrix_B[0][k][1]); n = ceil(Part_Matrix_B[0][k][1]);
		
		Area1 = (m - Part_Matrix_A[0][k][1]); Area2 = (Part_Matrix_A[0][k][1] - i);
		Area3 = (n - Part_Matrix_B[0][k][1]); Area4 = (Part_Matrix_B[0][k][1] - j);
		
		chargenodeA[0][i] += Area1;
		chargenodeA[0][m] += Area2;
		chargenodeB[0][j] += Area3;
		chargenodeB[0][n] += Area4;}
	
	chargenodeA[0][0] = chargenodeA[0][0] + chargenodeA[0][Lx];
	chargenodeA[0][Lx] = chargenodeA[0][0];
	chargenodeB[0][0] = chargenodeB[0][0] + chargenodeB[0][Lx];
	chargenodeB[0][Lx] = chargenodeB[0][0];
	
	double pcn = .5/(double)ppc;
	for (ic=0; ic<=Lx; ic++) {
		normchargenodeA[0][ic] = pcn*chargenodeA[0][ic];
		normchargenodeB[0][ic] = pcn*chargenodeB[0][ic];
		normdiff[0][ic] = 1 - normchargenodeA[0][ic] - normchargenodeB[0][ic];
		rho[ic] = normdiff[0][ic];}
	
	printf("\nTimestep:\t%4d\n",0);
	
	update_phi_field(rho);
	for (ef=0; ef<Lx; ef++) {
		phi[0][ef] = rho[ef];
		pef[ef] = phi[0][ef];}
	
	update_E_field(pef);
	for (ef=0; ef<Lx; ef++) {
		efield[0][ef] = pef[ef];}
	
	phi[0][Lx] = phi[0][0];
	efield[0][Lx] = efield[0][0];
	
}

// Update the particle position and velocity
// This is the part we need to parallelize
// 'periodic' means the periodic boundary condition for particles - if it's going to right-bound, it shows up in the left, vice versa
void periodic_move_node() {
    //Moved update node charge loop into part of the main loop.
    //parallelized with pthreads, required synchronization of chargeNodeB and chargenodeA variables via mutex locks
    //increased locality by using register variables
	int ii,kk,i,j,m,n,ic,ef;
	int ddx = 1/dx;
	double ef_inta,ef_intb;
	double xposa,xposb,xposnewa,xposnewb;
	double v_xa,v_xnewa,v_xb,v_xnewb,v_ya,v_za,v_yb,v_zb;
	double Area1 = 0, Area2 = 0, Area3 = 0, Area4 = 0;
	double pcn = .5/(double)ppc;
	double rho[Lx+1],pef[Lx+1];
    pthread_attr_t attr;
    pthread_t threads[NUM_THREADS];
	for (ii=1; ii<tt; ii++) {// For each time step - we can't parallelize the time-serial work...
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        int o = 0, l = 0;
        for (o = 0; o<NUM_THREADS; o++) {
            struct thread_data *thread_d = (struct thread_data *)malloc(sizeof(struct thread_data));
            thread_d->threadNum = o;
            thread_d->numElementsPerThread = Npart/NUM_THREADS;
            thread_d->numElements = Npart/NUM_THREADS;
            thread_d->currentii = ii;
            if (o == NUM_THREADS - 1) {
                thread_d->numElements += Npart - NUM_THREADS*Npart/NUM_THREADS;
            }
            int rc = pthread_create(&threads[o], &attr, periodic_move_parallel, (void *)thread_d);
            if (rc) {
                printf("unable to create thread");
            }
        }
        pthread_attr_destroy(&attr);
        for (l=0; l<NUM_THREADS; l++) {
            pthread_join(threads[l], NULL);
        }
        
		//	Periodic Node Charge
		chargenodeA[ii][0] = chargenodeA[ii][0] + chargenodeA[ii][Lx];
		chargenodeA[ii][Lx] = chargenodeA[ii][0];
		chargenodeB[ii][0] = chargenodeB[ii][0] + chargenodeB[ii][Lx];
		chargenodeB[ii][Lx] = chargenodeB[ii][0];
		
		//	Normalize Node Charge
		for (ic=0; ic<=Lx; ic++) {
			normchargenodeA[ii][ic] = pcn*chargenodeA[ii][ic];
			normchargenodeB[ii][ic] = pcn*chargenodeB[ii][ic];
			normdiff[ii][ic] = 1 - normchargenodeA[ii][ic] - normchargenodeB[ii][ic];
			rho[ic] = normdiff[ii][ic];
		}
		
		if (ii==.2*tt-1||ii ==.4*tt-1||ii ==.6*tt-1||ii==.8*tt-1||ii==tt-1) {printf("Timestep:\t%4d\n",ii + 1);}
		
		update_phi_field(rho);
		for (ef=0; ef<Lx; ef++) {
			phi[ii][ef] = rho[ef];
			pef[ef] = phi[ii][ef];
		}
		
		update_E_field(pef);
		for (ef=0; ef<Lx; ef++) {
			efield[ii][ef] = pef[ef];
		}
		
		phi[ii][Lx] = phi[ii][0];
		efield[ii][Lx] = efield[ii][0];
	}
	
}

void * periodic_move_parallel(void *thread_d) {
    int ddx = 1/dx;
    int kk = 0, ii = 0, i, j, m, n;
    double xposa,xposb,xposnewa,xposnewb, ef_inta, ef_intb, ef_save;
	double v_xa,v_xnewa,v_xb,v_xnewb,v_ya,v_za,v_yb,v_zb;
    double Area1 = 0, Area2 = 0, Area3 = 0, Area4 = 0;
    struct thread_data *mydata = (struct thread_data *) thread_d;
    int mystartindex = mydata->threadNum * mydata->numElementsPerThread;
    ii = mydata->currentii;
    for (kk=mystartindex; kk<mystartindex + mydata->numElements; kk++) { // @@@@@@@@@@@@@@@@@ Note the starting and ending index, 0 --> Npart, parallelize this
		//	Electron Population A
        xposa = Part_Matrix_A[ii-1][kk][1];
        v_xa = Part_Matrix_A[ii-1][kk][4];
        v_ya = Part_Matrix_A[ii-1][kk][5];
        v_za = Part_Matrix_A[ii-1][kk][6];
        
		//	Move Particles
        i = floor(xposa); m = ceil(xposa);
        ef_save = efield[ii-1][i];
        ef_inta = (efield[ii-1][m] - ef_save)*(xposa - i)*ddx + ef_save;
        v_xnewa = v_xa - ef_inta*dt;
        xposnewa = xposa + v_xnewa*dt;
        
		//	Update Particle Positions and Velocities
        Part_Matrix_A[ii][kk][0] = Part_Matrix_A[ii-1][kk][0];
        Part_Matrix_A[ii][kk][1] = xposnewa;
        Part_Matrix_A[ii][kk][4] = v_xnewa;
        Part_Matrix_A[ii][kk][5] = v_ya;
        Part_Matrix_A[ii][kk][6] = v_za;
        
		//	Periodic Particle Position
        if (xposnewa < 0) {
            xposnewa = xposnewa + Lx;
            Part_Matrix_A[ii][kk][1] = xposnewa;}
        if (xposnewa > Lx) {
            xposnewa = xposnewa - Lx;
            Part_Matrix_A[ii][kk][1] = xposnewa;}
		
		//	Electron Population B
        xposb = Part_Matrix_B[ii-1][kk][1];
        v_xb = Part_Matrix_B[ii-1][kk][4];
        v_yb = Part_Matrix_B[ii-1][kk][5];
        v_zb = Part_Matrix_B[ii-1][kk][6];
        
		//	Move Particles
        j = floor(xposb); n = ceil(xposb);
        ef_intb = (efield[ii-1][n] - efield[ii-1][j])*(xposb - j)*ddx + efield[ii-1][j];
        v_xnewb = v_xb - ef_intb*dt;
        xposnewb = xposb + v_xnewb*dt;
        
		//	Update Particle Positions and Velocities
        Part_Matrix_B[ii][kk][0] = Part_Matrix_B[ii-1][kk][0];
        Part_Matrix_B[ii][kk][1] = xposnewb;
        Part_Matrix_B[ii][kk][4] = v_xnewb;
        Part_Matrix_B[ii][kk][5] = v_yb;
        Part_Matrix_B[ii][kk][6] = v_zb;
        
		//	Periodic Particle Position
        if (xposnewb < 0) {
            xposnewb = xposnewb + Lx;
            Part_Matrix_B[ii][kk][1] = xposnewb;
        }
        if (xposnewb > Lx) {
            xposnewb = xposnewb - Lx;
            Part_Matrix_B[ii][kk][1] = xposnewb;
        }
        
        i = floor(xposnewa); m = ceil(xposnewa);
        j = floor(xposnewb); n = ceil(xposnewb);
        Area1 = (m - xposnewa); Area2 = (xposnewa - i);
        Area3 = (n - xposnewb); Area4 = (xposnewb - j);
        pthread_mutex_lock(&mutexcharge);
        chargenodeA[ii][i] += Area1;
        chargenodeA[ii][m] += Area2;
        chargenodeB[ii][j] += Area3;
        chargenodeB[ii][n] += Area4;
        pthread_mutex_unlock(&mutexcharge);

    }
    pthread_exit(NULL);
}

double update_phi_field(double rho_phi[]) {
	
	//	Construct B Matrix
	int i; double Bmatrix[mad];
	for (i=0; i<Lx; i++) {
		Bmatrix[i] = -delta2*rho_phi[i];}
	
	//	Implement SOR Solution for Phi
	double phi_new[NMAX] = {0}, phi_old[NMAX] = {0}, Matrix_Product[NMAX] = {0};
	double tol = 1e-6; double resid[NMAX] = {0}, resid_norm = 1;
	
	int k,l; double sum_mp, sum_norm;
	int counter = 0;
	while (resid_norm > tol && counter < 1000) {
		for (k=0; k<Lx; k++) {phi_old[k] = phi_new[k];}
		for (k=0; k<Lx; k++) {
			sum_mp = 0; sum_norm = 0;
			for (l=0; l<Lx; l++) {sum_mp += (double)Amatrix[k][l]*phi_new[l];}
			Matrix_Product[k] = sum_mp;}
		for (k=0; k<Lx; k++) {
			resid[k] = Matrix_Product[k] - Bmatrix[k];
			phi_new[k] = phi_old[k] + .5*omega*resid[k];
			sum_norm += resid[k]*resid[k];}
		resid_norm = sqrt(sum_norm);
		counter++;
		//	printf("Residuals:\t%f\n",resid_norm);
	}

	
	for (i=0; i<Lx; i++) {rho_phi[i] = phi_new[i];}
	
	return(1);
}

double update_E_field(double phi_ef[]) {

	//	Differentiate for E Field
	int i; double eff[Lx+1];
	double ddx = 1/dx;
	for (i=0; i<Lx; i++) {
		if (i == 0) {eff[i] = .5*ddx*(phi_ef[Lx-1] - phi_ef[1]);}
		else if (i == Lx-1) {eff[i] = .5*ddx*(phi_ef[Lx-2] - phi_ef[0]);}
		else {eff[i] = .5*ddx*(phi_ef[i-1] - phi_ef[i+1]);}}
	
	for (i=0; i<Lx; i++) {phi_ef[i] = eff[i];}
	
	return(1);
}

void output() {
	
	int usa = 125; int nn; double Bin[usa+1];
	double binsize = 50*v_th/usa; double norm = 1/(2*binsize*Npart);
	for (nn=0; nn<=usa; nn++){Bin[nn] = -25*v_th + binsize*nn;}
	
	int i,kk,wi;
	double histxia[mad]={0}, histyia[mad]={0}, histzia[mad]={0};
	double histxfa[mad]={0}, histyfa[mad]={0}, histzfa[mad]={0};
	double nhistxia[mad]={0}, nhistyia[mad]={0}, nhistzia[mad]={0};
	double nhistxfa[mad]={0}, nhistyfa[mad]={0}, nhistzfa[mad]={0};
	double histxib[mad]={0}, histyib[mad]={0}, histzib[mad]={0};
	double histxfb[mad]={0}, histyfb[mad]={0}, histzfb[mad]={0};
	double nhistxib[mad]={0}, nhistyib[mad]={0}, nhistzib[mad]={0};
	double nhistxfb[mad]={0}, nhistyfb[mad]={0}, nhistzfb[mad]={0};
	
	for (kk = 0; kk<Npart; kk++){
		for (wi=0; wi<usa; wi++) {
			if ((Part_Matrix_A[0][kk][4]>Bin[wi] - .5*binsize) && (Part_Matrix_A[0][kk][4]<Bin[wi+1] - .5*binsize)){
				histxia[wi]++;}
			if ((Part_Matrix_A[0][kk][5]>Bin[wi] - .5*binsize) && (Part_Matrix_A[0][kk][5]<Bin[wi+1] - .5*binsize)){
				histyia[wi]++;}
			if ((Part_Matrix_A[0][kk][6]>Bin[wi] - .5*binsize) && (Part_Matrix_A[0][kk][6]<Bin[wi+1] - .5*binsize)){
				histzia[wi]++;}
			if ((Part_Matrix_A[tt-1][kk][4]>Bin[wi] - .5*binsize) && (Part_Matrix_A[tt-1][kk][4]<Bin[wi+1] - .5*binsize)){
				histxfa[wi]++;}
			if ((Part_Matrix_A[tt-1][kk][5]>Bin[wi] - .5*binsize) && (Part_Matrix_A[tt-1][kk][5]<Bin[wi+1] - .5*binsize)){
				histyfa[wi]++;}
			if ((Part_Matrix_A[tt-1][kk][6]>Bin[wi] - .5*binsize) && (Part_Matrix_A[tt-1][kk][6]<Bin[wi+1] - .5*binsize)){
				histzfa[wi]++;}
			if ((Part_Matrix_B[0][kk][4]>Bin[wi] - .5*binsize) && (Part_Matrix_B[0][kk][4]<Bin[wi+1] - .5*binsize)){
				histxib[wi]++;}
			if ((Part_Matrix_B[0][kk][5]>Bin[wi] - .5*binsize) && (Part_Matrix_B[0][kk][5]<Bin[wi+1] - .5*binsize)){
				histyib[wi]++;}
			if ((Part_Matrix_B[0][kk][6]>Bin[wi] - .5*binsize) && (Part_Matrix_B[0][kk][6]<Bin[wi+1] - .5*binsize)){
				histzib[wi]++;}
			if ((Part_Matrix_B[tt-1][kk][4]>Bin[wi] - .5*binsize) && (Part_Matrix_B[tt-1][kk][4]<Bin[wi+1] - .5*binsize)){
				histxfb[wi]++;}
			if ((Part_Matrix_B[tt-1][kk][5]>Bin[wi] - .5*binsize) && (Part_Matrix_B[tt-1][kk][5]<Bin[wi+1] - .5*binsize)){
				histyfb[wi]++;}
			if ((Part_Matrix_B[tt-1][kk][6]>Bin[wi] - .5*binsize) && (Part_Matrix_B[tt-1][kk][6]<Bin[wi+1] - .5*binsize)){
				histzfb[wi]++;}}
	}
	
	for (i=0; i<=usa; i++){
		nhistxia[i] = norm*histxia[i]; nhistyia[i] = norm*histyia[i]; nhistzia[i] = norm*histzia[i];
		nhistxfa[i] = norm*histxfa[i]; nhistyfa[i] = norm*histyfa[i]; nhistzfa[i] = norm*histzfa[i];
		nhistxib[i] = norm*histxib[i]; nhistyib[i] = norm*histyib[i]; nhistzib[i] = norm*histzib[i];
		nhistxfb[i] = norm*histxfb[i]; nhistyfb[i] = norm*histyfb[i]; nhistzfb[i] = norm*histzfb[i];
	}
	
	FILE *file1;
	file1 = fopen("part_v_hist_initial.txt", "w");
	for (wi=0; wi<=usa; wi++){
		fprintf(file1,"%6.2f\t\t%f\t%f\t%f\t%f\t%f\t%f\n",Bin[wi],
				nhistxia[wi],nhistyia[wi],nhistzia[wi],nhistxib[wi],nhistyib[wi],nhistzib[wi]);}
	fclose(file1);
	
	FILE *file2;
	file2 = fopen("num_den_profile_initial.txt", "w");
	for (i=0; i<=Lx; i++) {
		fprintf(file2,"%d\t%9.6f\t%8.6f\t%8.6f\t%9.6f\n",i,
				phi[0][i],normchargenodeA[0][i],normchargenodeB[0][i],normdiff[0][i]);}
	fclose(file2);
	
	FILE *file3;
	file3 = fopen("part_a_xpos.txt", "w");
	for (kk=0; kk<Npart; kk++) {
		fprintf(file3,"%5.0f",Part_Matrix_A[0][kk][0]);
		for (i=0; i<tt; i+=tt/20) {
			fprintf(file3,"\t%6.3f",Part_Matrix_A[i][kk][1]);}
		fprintf(file3,"\n");}
	fclose(file3);
	
	FILE *file4;
	file4 = fopen("part_b_xpos.txt", "w");
	for (kk=0; kk<Npart; kk++) {
		fprintf(file4,"%5.0f",Part_Matrix_B[0][kk][0]);
		for (i=0; i<tt; i+=tt/20) {
			fprintf(file4,"\t%6.3f",Part_Matrix_B[i][kk][1]);}
		fprintf(file4,"\n");}
	fclose(file4);
	
	FILE *file5;
	file5 = fopen("part_v_hist_final.txt", "w");
	for (wi=0; wi<=usa; wi++){
		fprintf(file5,"%6.2f\t\t%f\t%f\t%f\t%f\t%f\t%f\n",Bin[wi],
				nhistxfa[wi],nhistyfa[wi],nhistzfa[wi],nhistxfb[wi],nhistyfb[wi],nhistzfb[wi]);}
	fclose(file5);
	
	FILE *file6;
	file6 = fopen("num_den_profile_final.txt", "w");
	for (i=0; i<=Lx; i++) {
		fprintf(file6,"%d\t%9.6f\t%8.6f\t%8.6f\t%9.6f\n",i,
				phi[tt-1][i],normchargenodeA[tt-1][i],normchargenodeB[tt-1][i],normdiff[tt-1][i]);}
	fclose(file6);
	
	FILE *file7;
	file7 = fopen("part_a_vx.txt", "w");
	for (kk=0; kk<Npart; kk++) {
		fprintf(file7,"%5.0f",Part_Matrix_A[0][kk][0]);
		for (i=0; i<tt; i+=tt/20) {
			fprintf(file7,"\t%7.3f",Part_Matrix_A[i][kk][4]);}
		fprintf(file7,"\n");}
	fclose(file7);
	
	FILE *file8;
	file8 = fopen("part_b_vx.txt", "w");
	for (kk=0; kk<Npart; kk++) {
		fprintf(file8,"%5.0f",Part_Matrix_B[0][kk][0]);
		for (i=0; i<tt; i+=tt/20) {
			fprintf(file8,"\t%7.3f",Part_Matrix_B[i][kk][4]);}
		fprintf(file8,"\n");}
	fclose(file8);
	
}

void vel_temp() {
	
	int ii,kk;
	double sum_vx, sqsum_vx;
	double avg_vxa[tt],dv2_xa[tt];
	double avg_vxb[tt],dv2_xb[tt];
	double Tempa[tt],Tempb[tt];
	for (ii=0;ii<tt;ii++) {
		sum_vx = 0; sqsum_vx = 0;
		for (kk=0; kk<Npart; kk++) {
			sum_vx += Part_Matrix_A[ii][kk][4];
			sqsum_vx += Part_Matrix_A[ii][kk][4]*Part_Matrix_A[ii][kk][4];}
		
		avg_vxa[ii] = sum_vx/Npart;
		dv2_xa[ii] = sqsum_vx/Npart - avg_vxa[ii]*avg_vxa[ii];
		
		sum_vx = 0; sqsum_vx = 0;
		for (kk=0; kk<Npart; kk++) {
			sum_vx += Part_Matrix_B[ii][kk][4];
			sqsum_vx += Part_Matrix_B[ii][kk][4]*Part_Matrix_B[ii][kk][4];}
		
		avg_vxb[ii] = sum_vx/Npart;
		dv2_xb[ii] = sqsum_vx/Npart - avg_vxb[ii]*avg_vxb[ii];
		
		Tempa[ii] = 2*dv2_xa[ii]; Tempb[ii] = 2*dv2_xb[ii];
		
		if (ii == 0){
			printf("\n\n=========Initial=========\n\n");
			printf("avg_vx:\t%6.3f\t%8.3f\n",avg_vxa[ii],avg_vxb[ii]);
			printf("dv2_x:\t\t%6.3f\t%8.3f\n",dv2_xa[ii],dv2_xb[ii]);
			printf("\nTemps:\t\t%6.3f\t%8.3f",Tempa[ii],Tempb[ii]);}
		
		if (ii == tt-1){
			printf("\n\n==========Final==========\n\n");
			printf("avg_vx:\t%6.3f\t%8.3f\n",avg_vxa[ii],avg_vxb[ii]);
			printf("dv2_x:\t\t%6.3f\t%8.3f\n",dv2_xa[ii],dv2_xb[ii]);
			printf("\nTemps:\t\t%6.3f\t%8.3f",Tempa[ii],Tempb[ii]);
			printf("\n\n=========================\n\n");}
	}
	
	FILE *file1;
	file1 = fopen("velocity_temp_hist.txt", "w");
	for (ii=0; ii<tt; ii++) {
		fprintf(file1,"%4d\t%f\t%f\t%f\t%f\n",ii,
				avg_vxa[ii],Tempa[ii],avg_vxb[ii],Tempb[ii]);}
	fclose(file1);
	
}

void energy() {
	
	int i,k,ef;
	double v_x2[2][tt] = {0}, v_y2[2][tt] = {0},v_z2[2][tt] = {0}, ef2[tt] = {0};
	double v_tot2 [2][tt] = {0};
	double ap = (double) .5*Lx/(ppc*Lx);
	
	for (i=0; i<tt; i++) {
		for (k=0; k<Npart; k++) {
			v_x2[0][i] += Part_Matrix_A[i][k][4]*Part_Matrix_A[i][k][4];
			v_y2[0][i] += Part_Matrix_A[i][k][5]*Part_Matrix_A[i][k][5];
			v_z2[0][i] += Part_Matrix_A[i][k][6]*Part_Matrix_A[i][k][6];
			v_tot2[0][i] = .5*(v_x2[0][i] + v_y2[0][i] + v_y2[0][i]);
			v_x2[1][i] += Part_Matrix_B[i][k][4]*Part_Matrix_B[i][k][4];
			v_y2[1][i] += Part_Matrix_B[i][k][5]*Part_Matrix_B[i][k][5];
			v_z2[1][i] += Part_Matrix_B[i][k][6]*Part_Matrix_B[i][k][6];
			v_tot2[1][i] = .5*(v_x2[1][i] + v_y2[1][i] + v_y2[1][i]);}
		
		for (ef=0; ef<=Lx; ef++) {
			ef2[i] += .5*efield[i][ef]*efield[i][ef];}
	}
	
	FILE *file1;
	file1 = fopen("energy.txt", "w");
	for (i=0; i<tt; i++) {
		fprintf(file1,"%3d\t%f\t%f\t%f\n",i,
				ap*v_tot2[0][i],ap*v_tot2[1][i],ef2[i]);}
	fclose(file1);
	
}

double rng(double num) {
	
	int MAXVALUE;
	double nMAXVALUE;
	MAXVALUE = 100000;
	nMAXVALUE = (double)1/MAXVALUE;
	int value;
	
	value = rand() % (MAXVALUE + 1);
	if (value > MAXVALUE)
		value = 0;
	
	num = value*nMAXVALUE;
	
	return(num);
}
