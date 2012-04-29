// Hybrid PIC Plasma Sheath

#include <stdio.h>
#include <stdlib.h>
#include <libc.h>
#include <math.h>
#include "ps.h"
#include <pthread.h>
double elapsed, seconds;
#define NUM_THREADS 1
pthread_mutex_t mutex_move = PTHREAD_MUTEX_INITIALIZER;
int lastk = 0;
struct thread_data {
    int threadNum;
    int numElements;
    int numElementsPerThread;
    int currentii;
    int numIterations;
};
int main () {
	
	struct timeval tv;
	struct timezone tz;
	struct tm*tm;
	tm = localtime(&tv.tv_sec);
	gettimeofday(&tv, &tz);
	//srand(tv.tv_usec * tv.tv_sec);
	
	v_o = 4; v_th = sqrt(.2);
	
	//Lx = 40; Ly = 20;
	Lx = 40; Ly = 20;
    dx = 1; dy = 1; delta2 = dx*dx;
	imax = Lx - 1; jmax = Ly;
	nmax = imax*jmax;
	ppc = 1200;
    //ppc = 200; 
    LeftBC = 0; RightBC = -10;
    pcn = 1/(double)ppc;
	Npart = Lx*Ly*ppc;
	
	pi = 4*atan(1);
	nt = .5;
	dt = nt/(v_o*(1 + erf(v_th/v_o)));
	//finaltime = 50;
	finaltime = 20;
    omega = 1;
	initial_conditions();
    seconds = read_timer();
	inject_node_move();
    elapsed = read_timer() - seconds;
	outputs();
	
	FILE *file1;
	file1 = fopen("Sys_Cond.txt", "w");
	fprintf(file1, "%d\t%d\t%d\t%d\t%d",Lx,Ly,dx,dy,finaltime);
	fclose(file1);
	printf("Elapsed time (parallel): %f seconds\n", elapsed);
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

void initial_conditions() {
	
	int i,j,k;
	double ni_phi[AMIN];
	
	//	Initiate Ion Node Charge & Normalize
	for (i=0; i<=Lx; i++) {
		for (j=0; j<=Ly; j++) {
			ion_node[i][j] = 0;
			norm_inode[i][j] = pcn*ion_node[i][j];}
	}
	
	//	Disperse Ion Node into a Single Array
	for (i=0; i<imax; i++) {
		for (j=0; j<jmax; j++) {
			k = imax*j + i;
			ni_phi[k] = norm_inode[i+1][j];}
	}
	
	//	Construct Amatrix
	for (i=0; i<imax; i++) {
		for (j=0; j<jmax; j++) {
			k = imax*j + i;
			if (i == 0 && j == 0) {
				Amatrix[k][k+1] = 1;
				Amatrix[k][k+imax] = 1;
				Amatrix[k][k+(jmax-1)*imax] = 1;}
			else if (i == imax-1 && j == 0) {
				Amatrix[k][k-1] = 1;
				Amatrix[k][k+imax] = 1;
				Amatrix[k][k+(jmax-1)*imax] = 1;}
			else if (i != 0 && i != imax-1 && j == 0) {
				Amatrix[k][k-1] = 1; Amatrix[k][k+1] = 1;
				Amatrix[k][k+imax] = 1;
				Amatrix[k][k+(jmax-1)*imax] = 1;}
			else if (i == 0 && j == jmax-1) {
				Amatrix[k][k+1] = 1;
				Amatrix[k][k-imax] = 1;
				Amatrix[k][k-(jmax-1)*imax] = 1;}
			else if (i == imax-1 && j == jmax-1) {
				Amatrix[k][k-1] = 1;
				Amatrix[k][k-imax] = 1;
				Amatrix[k][k-(jmax-1)*imax] = 1;}
			else if (i != 0 && i != imax-1 && j == jmax-1) {
				Amatrix[k][k-1] = 1; Amatrix[k][k+1] = 1;
				Amatrix[k][k-imax] = 1;
				Amatrix[k][k-(jmax-1)*imax] = 1;}
			else if (i == 0 && j != 0) {
				Amatrix[k][k+1] = 1;
				Amatrix[k][k-imax] = 1; Amatrix[k][k+imax] = 1;}
			else if (i == imax-1 && j != jmax-1) {
				Amatrix[k][k-1] = 1;
				Amatrix[k][k-imax] = 1; Amatrix[k][k+imax] = 1;}
			else {
				Amatrix[k][k-1] = 1; Amatrix[k][k+1] = 1;
				Amatrix[k][k-imax] = 1; Amatrix[k][k+imax] = 1;}}
	}
	
	printf("\nSys Length:%5d  x %3d",Lx,Ly);
	printf("\nCell Size:%7.1f x %3.1f",(double) dx,(double) dy);
	printf("\nAmatrix:%9d x %3d\n",nmax,nmax);
	
	printf("\n=======================");
	printf("\n    Time Step:%5d",0);
	printf("\n=======================\n");
	
	//	Solve and Map Phi
	nl_phi_solver(ni_phi);
	for (i=0; i<imax; i++) {
		for (j=0; j<jmax; j++) {
			k = imax*j + i;
			phi[0][j] = LeftBC;
			phi[i+1][j] = ni_phi[k];
			phi[imax+1][j] = RightBC;}
	}
	for (i=0; i<=Lx; i++) {
		phi[i][Ly] = phi[i][0];
	}
	
	//	Differentiate to Obtain E Field
	for (j=0; j<=Ly; j++) {
		for (i=1; i<Lx; i++) {
			Ex[i][j] = -.5*(phi[i+1][j] - phi[i-1][j])/dx;}
		Ex[0][j] = 2*(phi[0][j] - phi[1][j])/dx - Ex[1][j];
		Ex[Lx][j] = 2*(phi[Lx-1][j] - phi[Lx][j])/dx - Ex[Lx-1][j];
	}
	for (i=0; i<=Lx; i++) {
		for (j=1; j<Ly; j++) {
			Ey[i][j] = -.5*(phi[i][j+1] - phi[i][j-1])/dy;}
		Ey[i][0] = -.5*(phi[i][1] - phi[i][Ly-1])/dy;
		Ey[i][Ly] = Ey[i][0];
	}
	
	FILE *file1;
	file1 = fopen("PS_Phi_Initial.txt", "w");
	for (j=0; j<=Ly; j++) {
		for (i=0; i<=Lx; i++) {
			fprintf(file1, "%10.6f\t",phi[i][j]);}
		fprintf(file1, "\n");}
	fclose(file1);
	
}

void inject_node_move() {
	
	int ii,inj,apd,idx;
	int injcount;
	int i,j,k,m,n;
	int Pinj = nt*Ly*ppc;
	int particletracker =  0;
	double a,vrand1,vrand2,th1,th2;
	double x_pos,y_pos,z_pos;
	double v_x,v_y,v_z;
	double Area1,Area2,Area3,Area4;
	double Ex1,Ex2,Ey1,Ey2,Ex_int,Ey_int;
	double xposnew,yposnew;
	double vxnew,vynew;
	double ni_phi[AMIN];
    pthread_attr_t attr;
    pthread_t threads[NUM_THREADS];
	for (ii=1; ii<=finaltime; ii++) {
        
        //Regular loop: No optimization
		/*for (k=0; k<newcount; k++) {
			for (idx=0; idx<indx; idx++) {
				PMO[k][idx] = Part_Matrix[k][idx];}
		}*/
        
        //Tried loop tiling of the 2-D loop above with many, many different tile sizes, no improvement
        /*int kk=0;
        int iidx = 0;
        int tile1 = 128;
        int tile2 = 7;
        for (k=0; k<newcount; k+=tile1) {
			for (idx=0; idx<indx; idx+=tile2) {
                for (kk=k; kk<(((k + tile1) < newcount)  ? k+tile1 : newcount); kk++) {
                    for (iidx = idx; iidx<(((idx + tile2) < indx) ? idx + tile2: indx); iidx++) {
                        PMO[kk][iidx] = Part_Matrix[kk][iidx];
                    }
                }
            }
		}*/
        
        //Settled for unrolling the entirety of the inner loop and 2 iterations in the outer one. Gave a little speedup (~1 second)
        for (k=0; k<newcount; k+=2) {
            PMO[k][0] = Part_Matrix[k][0];
            PMO[k][1] = Part_Matrix[k][1];
            PMO[k][2] = Part_Matrix[k][2];
            PMO[k][3] = Part_Matrix[k][3];
            PMO[k][4] = Part_Matrix[k][4];
            PMO[k][5] = Part_Matrix[k][5];
            PMO[k][6] = Part_Matrix[k][6];
            if (k + 1 < newcount) {
                PMO[k+1][0] = Part_Matrix[k+1][0];
                PMO[k+1][1] = Part_Matrix[k+1][1];
                PMO[k+1][2] = Part_Matrix[k+1][2];
                PMO[k+1][3] = Part_Matrix[k+1][3];
                PMO[k+1][4] = Part_Matrix[k+1][4];
                PMO[k+1][5] = Part_Matrix[k+1][5];
                PMO[k+1][6] = Part_Matrix[k+1][6];
                /*if (k +2 < newcount) {
                    PMO[k+2][0] = Part_Matrix[k+2][0];
                    PMO[k+2][1] = Part_Matrix[k+2][1];
                    PMO[k+2][2] = Part_Matrix[k+2][2];
                    PMO[k+2][3] = Part_Matrix[k+2][3];
                    PMO[k+2][4] = Part_Matrix[k+2][4];
                    PMO[k+2][5] = Part_Matrix[k+2][5];
                    PMO[k+2][6] = Part_Matrix[k+2][6];
                    if (k+3 < newcount) {
                        PMO[k+3][0] = Part_Matrix[k+3][0];
                        PMO[k+3][1] = Part_Matrix[k+3][1];
                        PMO[k+3][2] = Part_Matrix[k+3][2];
                        PMO[k+3][3] = Part_Matrix[k+3][3];
                        PMO[k+3][4] = Part_Matrix[k+3][4];
                        PMO[k+3][5] = Part_Matrix[k+3][5];
                        PMO[k+3][6] = Part_Matrix[k+3][6];
                        if (k+4 < newcount) {
                            PMO[k+4][0] = Part_Matrix[k+4][0];
                            PMO[k+4][1] = Part_Matrix[k+4][1];
                            PMO[k+4][2] = Part_Matrix[k+4][2];
                            PMO[k+4][3] = Part_Matrix[k+4][3];
                            PMO[k+4][4] = Part_Matrix[k+4][4];
                            PMO[k+4][5] = Part_Matrix[k+4][5];
                            PMO[k+4][6] = Part_Matrix[k+4][6];
                        }
                    }
                }*/
            }
         }
        
		//	Inject Particles @ x = 0
		for (inj=0; inj<Pinj; inj+=10) {
			particletracker++;
			vrand1 = sqrt(-log(1 - .999999*rng(a)));
			vrand2 = sqrt(-log(1 - .999999*rng(a)));
			th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
			
			x_pos = nt*rng(a);
			y_pos = Ly*rng(a);
			z_pos = 0;
			v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
			v_y = v_th*vrand1*sin(th1);
			v_z = v_th*vrand2*cos(th2);
			
			//dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
            dpm[inj][0] = particletracker;
			dpm[inj][1] = x_pos;
			dpm[inj][2] = y_pos;
			dpm[inj][3] = z_pos;
			dpm[inj][4] = v_x;
			dpm[inj][5] = v_y;
			dpm[inj][6] = v_z;
            
            if (inj + 1 < Pinj) {
                particletracker++;
                vrand1 = sqrt(-log(1 - .999999*rng(a)));
                vrand2 = sqrt(-log(1 - .999999*rng(a)));
                th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                
                x_pos = nt*rng(a);
                y_pos = Ly*rng(a);
                z_pos = 0;
                v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                v_y = v_th*vrand1*sin(th1);
                v_z = v_th*vrand2*cos(th2);
                int inj_next = inj + 1;
                //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                dpm[inj_next][0] = particletracker;
                dpm[inj_next][1] = x_pos;
                dpm[inj_next][2] = y_pos;
                dpm[inj_next][3] = z_pos;
                dpm[inj_next][4] = v_x;
                dpm[inj_next][5] = v_y;
                dpm[inj_next][6] = v_z;
                
                if (inj + 2 < Pinj) {
                    particletracker++;
                    vrand1 = sqrt(-log(1 - .999999*rng(a)));
                    vrand2 = sqrt(-log(1 - .999999*rng(a)));
                    th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                    
                    x_pos = nt*rng(a);
                    y_pos = Ly*rng(a);
                    z_pos = 0;
                    v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                    v_y = v_th*vrand1*sin(th1);
                    v_z = v_th*vrand2*cos(th2);
                    int inj_next = inj + 2;
                    //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                    dpm[inj_next][0] = particletracker;
                    dpm[inj_next][1] = x_pos;
                    dpm[inj_next][2] = y_pos;
                    dpm[inj_next][3] = z_pos;
                    dpm[inj_next][4] = v_x;
                    dpm[inj_next][5] = v_y;
                    dpm[inj_next][6] = v_z;
                    if (inj + 3 < Pinj) {
                        particletracker++;
                        vrand1 = sqrt(-log(1 - .999999*rng(a)));
                        vrand2 = sqrt(-log(1 - .999999*rng(a)));
                        th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                        
                        x_pos = nt*rng(a);
                        y_pos = Ly*rng(a);
                        z_pos = 0;
                        v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                        v_y = v_th*vrand1*sin(th1);
                        v_z = v_th*vrand2*cos(th2);
                        int inj_next = inj + 3;
                        //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                        dpm[inj_next][0] = particletracker;
                        dpm[inj_next][1] = x_pos;
                        dpm[inj_next][2] = y_pos;
                        dpm[inj_next][3] = z_pos;
                        dpm[inj_next][4] = v_x;
                        dpm[inj_next][5] = v_y;
                        dpm[inj_next][6] = v_z;
                        if (inj + 4 < Pinj) {
                            particletracker++;
                            vrand1 = sqrt(-log(1 - .999999*rng(a)));
                            vrand2 = sqrt(-log(1 - .999999*rng(a)));
                            th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                            
                            x_pos = nt*rng(a);
                            y_pos = Ly*rng(a);
                            z_pos = 0;
                            v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                            v_y = v_th*vrand1*sin(th1);
                            v_z = v_th*vrand2*cos(th2);
                            int inj_next = inj + 4;
                            //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                            dpm[inj_next][0] = particletracker;
                            dpm[inj_next][1] = x_pos;
                            dpm[inj_next][2] = y_pos;
                            dpm[inj_next][3] = z_pos;
                            dpm[inj_next][4] = v_x;
                            dpm[inj_next][5] = v_y;
                            dpm[inj_next][6] = v_z;
                            if (inj + 5 < Pinj) {
                                particletracker++;
                                vrand1 = sqrt(-log(1 - .999999*rng(a)));
                                vrand2 = sqrt(-log(1 - .999999*rng(a)));
                                th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                                
                                x_pos = nt*rng(a);
                                y_pos = Ly*rng(a);
                                z_pos = 0;
                                v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                                v_y = v_th*vrand1*sin(th1);
                                v_z = v_th*vrand2*cos(th2);
                                int inj_next = inj + 5;
                                //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                                dpm[inj_next][0] = particletracker;
                                dpm[inj_next][1] = x_pos;
                                dpm[inj_next][2] = y_pos;
                                dpm[inj_next][3] = z_pos;
                                dpm[inj_next][4] = v_x;
                                dpm[inj_next][5] = v_y;
                                dpm[inj_next][6] = v_z;
                                if (inj + 6 < Pinj) {
                                    particletracker++;
                                    vrand1 = sqrt(-log(1 - .999999*rng(a)));
                                    vrand2 = sqrt(-log(1 - .999999*rng(a)));
                                    th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                                    
                                    x_pos = nt*rng(a);
                                    y_pos = Ly*rng(a);
                                    z_pos = 0;
                                    v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                                    v_y = v_th*vrand1*sin(th1);
                                    v_z = v_th*vrand2*cos(th2);
                                    int inj_next = inj + 6;
                                    //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                                    dpm[inj_next][0] = particletracker;
                                    dpm[inj_next][1] = x_pos;
                                    dpm[inj_next][2] = y_pos;
                                    dpm[inj_next][3] = z_pos;
                                    dpm[inj_next][4] = v_x;
                                    dpm[inj_next][5] = v_y;
                                    dpm[inj_next][6] = v_z;
                                    if (inj + 7 < Pinj) {
                                        particletracker++;
                                        vrand1 = sqrt(-log(1 - .999999*rng(a)));
                                        vrand2 = sqrt(-log(1 - .999999*rng(a)));
                                        th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                                        
                                        x_pos = nt*rng(a);
                                        y_pos = Ly*rng(a);
                                        z_pos = 0;
                                        v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                                        v_y = v_th*vrand1*sin(th1);
                                        v_z = v_th*vrand2*cos(th2);
                                        int inj_next = inj + 7;
                                        //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                                        dpm[inj_next][0] = particletracker;
                                        dpm[inj_next][1] = x_pos;
                                        dpm[inj_next][2] = y_pos;
                                        dpm[inj_next][3] = z_pos;
                                        dpm[inj_next][4] = v_x;
                                        dpm[inj_next][5] = v_y;
                                        dpm[inj_next][6] = v_z;
                                        if (inj + 8 < Pinj) {
                                            particletracker++;
                                            vrand1 = sqrt(-log(1 - .999999*rng(a)));
                                            vrand2 = sqrt(-log(1 - .999999*rng(a)));
                                            th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                                            
                                            x_pos = nt*rng(a);
                                            y_pos = Ly*rng(a);
                                            z_pos = 0;
                                            v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                                            v_y = v_th*vrand1*sin(th1);
                                            v_z = v_th*vrand2*cos(th2);
                                            int inj_next = inj + 8;
                                            //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                                            dpm[inj_next][0] = particletracker;
                                            dpm[inj_next][1] = x_pos;
                                            dpm[inj_next][2] = y_pos;
                                            dpm[inj_next][3] = z_pos;
                                            dpm[inj_next][4] = v_x;
                                            dpm[inj_next][5] = v_y;
                                            dpm[inj_next][6] = v_z;
                                            if (inj + 9 < Pinj) {
                                                particletracker++;
                                                vrand1 = sqrt(-log(1 - .999999*rng(a)));
                                                vrand2 = sqrt(-log(1 - .999999*rng(a)));
                                                th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
                                                
                                                x_pos = nt*rng(a);
                                                y_pos = Ly*rng(a);
                                                z_pos = 0;
                                                v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
                                                v_y = v_th*vrand1*sin(th1);
                                                v_z = v_th*vrand2*cos(th2);
                                                int inj_next = inj + 9;
                                                //dpm[inj][0] = (ii - 1)*Pinj+inj + 1;
                                                dpm[inj_next][0] = particletracker;
                                                dpm[inj_next][1] = x_pos;
                                                dpm[inj_next][2] = y_pos;
                                                dpm[inj_next][3] = z_pos;
                                                dpm[inj_next][4] = v_x;
                                                dpm[inj_next][5] = v_y;
                                                dpm[inj_next][6] = v_z;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
		}
        
        /*pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        int o = 0, l = 0;
        for (o = 0; o<NUM_THREADS; o++) {
            struct thread_data *thread_d = (struct thread_data *)malloc(sizeof(struct thread_data));
            thread_d->threadNum = o;
            thread_d->numElementsPerThread = Pinj/NUM_THREADS;
            thread_d->numElements = Pinj/NUM_THREADS;
            thread_d->currentii = ii;
            thread_d->numIterations = Pinj;
            if (o == NUM_THREADS - 1) {
                thread_d->numElements += Pinj - NUM_THREADS*(Pinj/NUM_THREADS);
            }
            int rc = pthread_create(&threads[o], &attr, inject_particles_parallel, (void *)thread_d);
            if (rc) {
                printf("unable to create thread");
            }
        }
        pthread_attr_destroy(&attr);
        for (l=0; l<NUM_THREADS; l++) {
            pthread_join(threads[l], NULL);
        }*/
		
        
		//	Append Injected Particles
		/*for (inj=0; inj<Pinj; inj++) {
			apd = newcount + inj;
			for (idx=0; idx<indx; idx++) {
				PMO[apd][idx] = dpm[inj][idx];}
		}*/
        
        for (inj=0; inj<Pinj; inj++) {
			 apd = newcount + inj;
             PMO[apd][0] = dpm[inj][0];
             PMO[apd][1] = dpm[inj][1];
             PMO[apd][2] = dpm[inj][2];
             PMO[apd][3] = dpm[inj][3];
             PMO[apd][4] = dpm[inj][4];
             PMO[apd][5] = dpm[inj][5];
             PMO[apd][6] = dpm[inj][6];
		}

		//	Update Ion Node Charge
		for (i=0; i<=Lx; i++) {
			for (j=0; j<=Ly; j++) {
				ion_node[i][j] = 0;}
		}
		injcount = newcount + Pinj;
		for (k=0; k<injcount; k++) {
			i = floor(PMO[k][1]);	m = ceil(PMO[k][1]);
			j = floor(PMO[k][2]);	n = ceil(PMO[k][2]);

			Area1 = (m - PMO[k][1])*(n - PMO[k][2]);
			Area2 = (PMO[k][1] - i)*(n - PMO[k][2]);
			Area3 = (PMO[k][1] - i)*(PMO[k][2] - j);
			Area4 = (m - PMO[k][1])*(PMO[k][2] - j);
			
			ion_node[i][j] += Area1;
			ion_node[m][j] += Area2;
			ion_node[m][n] += Area3;
			ion_node[i][n] += Area4;
		}
		//	Periodic in y
		for (i=0; i<=Lx; i++) {
			ion_node[i][0] += ion_node[i][Ly];
			ion_node[i][Ly] = ion_node[i][0];
		}
		
		//	Normalize Node Charge
		for (i=0; i<=Lx; i++) {
			for (j=0; j<=Ly; j++) {
				norm_inode[i][j] = pcn*ion_node[i][j];}
		}
		for (j=0; j<=Ly; j++) {
			norm_inode[0][j] = 2*pcn*ion_node[0][j];
			norm_inode[Lx][j] = 2*pcn*ion_node[Lx][j];
		}
		
		//	Disperse Ion Node into a Single Array
		for (i=0; i<imax; i++) {
			for (j=0; j<jmax; j++) {
				k = imax*j + i;
				ni_phi[k] = norm_inode[i+1][j];}
		}
		
		printf("\n=======================");
		printf("\n    Time Step:%5d",ii);
		printf("\n=======================\n");
		
		//	Solve and Map Phi
		nl_phi_solver(ni_phi);
		for (i=0; i<imax; i++) {
			for (j=0; j<jmax; j++) {
				k = imax*j + i;
				phi[0][j] = LeftBC;
				phi[i+1][j] = ni_phi[k];
				phi[imax+1][j] = RightBC;}
		}
		for (i=0; i<=Lx; i++) {
			phi[i][Ly] = phi[i][0];
		}
		
		//	Differentiate to Obtain E Field
		for (j=0; j<=Ly; j++) {
			for (i=1; i<Lx; i++) {
				Ex[i][j] = -.5*(phi[i+1][j] - phi[i-1][j])/dx;}
			Ex[0][j] = 2*(phi[0][j] - phi[1][j])/dx - Ex[1][j];
			Ex[Lx][j] = 2*(phi[Lx-1][j] - phi[Lx][j])/dx - Ex[Lx-1][j];
		}
		for (i=0; i<=Lx; i++) {
			for (j=1; j<Ly; j++) {
				Ey[i][j] = -.5*(phi[i][j+1] - phi[i][j-1])/dy;}
			Ey[i][0] = -.5*(phi[i][1] - phi[i][Ly-1])/dy;
			Ey[i][Ly] = Ey[i][0];
		}

		//	Move Particles
		for (k=0; k<injcount; k++) {
			x_pos = PMO[k][1];
			y_pos = PMO[k][2];
			z_pos = PMO[k][3];
			v_x = PMO[k][4];
			v_y = PMO[k][5];
			v_z = PMO[k][6];
			
			i = floor(x_pos); m = ceil(x_pos);
			j = floor(y_pos); n = ceil(y_pos);
			
			Ex1 = Ex[i][j] + (x_pos - i)*(Ex[m][j] - Ex[i][j])/dx;
			Ex2 = Ex[i][n] + (x_pos - i)*(Ex[m][n] - Ex[i][n])/dx;
			Ex_int = Ex1 + (y_pos - j)*(Ex2 - Ex1)/dy;
			// xposnew = x_pos + v_x*dt + .5*Ex_int*dt*dt;
			vxnew = v_x + Ex_int*dt;
			xposnew = x_pos + vxnew*dt;
			
			Ey1 = Ey[i][j] + (x_pos - i)*(Ey[m][j] - Ey[i][j])/dx;
			Ey2 = Ey[i][n] + (x_pos - i)*(Ey[m][n] - Ey[i][n])/dx;
			Ey_int = Ey1 + (y_pos - j)*(Ey2 - Ey1)/dy;
			// yposnew = y_pos + v_y*dt + .5*Ey_int*dt*dt;
			vynew = v_y + Ey_int*dt;
			yposnew = y_pos + vynew*dt;
			
			//	period boundary in y direction
			if (yposnew < 0) {
				yposnew += Ly;}
			if (yposnew > Ly) {
				yposnew += -Ly;}
			//	absorption boundary @ x = Lx
			if (xposnew > 0 && xposnew < Lx) {
				dpm[pct[ii]][0] = PMO[k][0];
				dpm[pct[ii]][1] = xposnew;
				dpm[pct[ii]][2] = yposnew;
				dpm[pct[ii]][3] = z_pos;
				dpm[pct[ii]][4] = vxnew;
				dpm[pct[ii]][5] = vynew;
				dpm[pct[ii]][6] = v_z;
				pct[ii]++;}
		}
        
        
        /*pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        int o = 0, l = 0;
        for (o = 0; o<NUM_THREADS; o++) {
            struct thread_data *thread_d = (struct thread_data *)malloc(sizeof(struct thread_data));
            thread_d->threadNum = o;
            thread_d->numElementsPerThread = injcount/NUM_THREADS;
            thread_d->numElements = injcount/NUM_THREADS;
            thread_d->currentii = ii;
            thread_d->numIterations = injcount;
            if (o == NUM_THREADS - 1) {
                thread_d->numElements += injcount - NUM_THREADS*(injcount/NUM_THREADS);
            }
            int rc = pthread_create(&threads[o], &attr, inject_node_move_parallel, (void *)thread_d);
            if (rc) {
                printf("unable to create thread");
            }
        }
        pthread_attr_destroy(&attr);
        for (l=0; l<NUM_THREADS; l++) {
            pthread_join(threads[l], NULL);
        }*/
        
		newcount = pct[ii];
		//	Update Particle Matrix
		for (k=0; k<newcount; k++) {
			for (idx=0; idx<indx; idx++) {
				Part_Matrix[k][idx] = dpm[k][idx];}
		}
		
	}
	
	printf("\nParticles Injected:\t%8.0f\n",Part_Matrix[k-1][0]);
	printf("Remaining Particles:\t%8d\n",newcount);

}

void * inject_node_move_parallel(void *thread_d) {
    int ii;
	int i,j,k,m,n;
	double vrand1,vrand2,th1,th2;
	double x_pos,y_pos,z_pos;
	double v_x,v_y,v_z;
	double Ex1,Ex2,Ey1,Ey2,Ex_int,Ey_int;
	double xposnew,yposnew;
	double vxnew,vynew;
    struct thread_data *mydata = (struct thread_data *) thread_d;
    int mystartindex = mydata->threadNum * mydata->numElementsPerThread;
    ii = mydata->currentii;
    for (k=mystartindex; k<mystartindex + mydata->numElements; k++) {
			x_pos = PMO[k][1];
			y_pos = PMO[k][2];
			z_pos = PMO[k][3];
			v_x = PMO[k][4];
			v_y = PMO[k][5];
			v_z = PMO[k][6];
			
			i = floor(x_pos); m = ceil(x_pos);
			j = floor(y_pos); n = ceil(y_pos);
			
			Ex1 = Ex[i][j] + (x_pos - i)*(Ex[m][j] - Ex[i][j])/dx;
			Ex2 = Ex[i][n] + (x_pos - i)*(Ex[m][n] - Ex[i][n])/dx;
			Ex_int = Ex1 + (y_pos - j)*(Ex2 - Ex1)/dy;
			// xposnew = x_pos + v_x*dt + .5*Ex_int*dt*dt;
			vxnew = v_x + Ex_int*dt;
			xposnew = x_pos + vxnew*dt;
			
			Ey1 = Ey[i][j] + (x_pos - i)*(Ey[m][j] - Ey[i][j])/dx;
			Ey2 = Ey[i][n] + (x_pos - i)*(Ey[m][n] - Ey[i][n])/dx;
			Ey_int = Ey1 + (y_pos - j)*(Ey2 - Ey1)/dy;
			// yposnew = y_pos + v_y*dt + .5*Ey_int*dt*dt;
			vynew = v_y + Ey_int*dt;
			yposnew = y_pos + vynew*dt;
			
			//	period boundary in y direction
			if (yposnew < 0) {
				yposnew += Ly;}
			if (yposnew > Ly) {
				yposnew += -Ly;}
			//	absorption boundary @ x = Lx
			if (xposnew > 0 && xposnew < Lx) {
                pthread_mutex_lock(&mutex_move);
                int s = pct[ii];
                dpm[s][0] = PMO[k][0];
                dpm[s][1] = xposnew;
                dpm[s][2] = yposnew;
                dpm[s][3] = z_pos;
                dpm[s][4] = vxnew;
                dpm[s][5] = vynew;
                dpm[s][6] = v_z;
                s++;
				pct[ii] = s;
                pthread_mutex_unlock(&mutex_move);
            
            }
        
    }
    pthread_exit(NULL);
}

void * inject_particles_parallel(void *thread_d) {
    int ii;
	int inj;
    double a;
	double vrand1,vrand2,th1,th2;
	double x_pos,y_pos,z_pos;
	double v_x,v_y,v_z;
    struct thread_data *mydata = (struct thread_data *) thread_d;
    int mystartindex = mydata->threadNum * mydata->numElementsPerThread;
    ii = mydata->currentii;
    for (inj=mystartindex; inj<mystartindex + mydata->numElements; inj++) {
     vrand1 = sqrt(-log(1 - .999999*rng(a)));
     vrand2 = sqrt(-log(1 - .999999*rng(a)));
     th1 = 2*pi*rng(a); th2 = 2*pi*rng(a);
     
     x_pos = nt*rng(a);
     y_pos = Ly*rng(a);
     z_pos = 0;
     v_x = v_o*(1 + erf(v_th/v_o)) + v_th*vrand1*cos(th1);
     v_y = v_th*vrand1*sin(th1);
     v_z = v_th*vrand2*cos(th2);
     dpm[inj][0] = (ii - 1)*mydata->numIterations+inj + 1;
     dpm[inj][1] = x_pos;
     dpm[inj][2] = y_pos;
     dpm[inj][3] = z_pos;
     dpm[inj][4] = v_x;
     dpm[inj][5] = v_y;
     dpm[inj][6] = v_z;
    }
    pthread_exit(NULL);
}


void outputs() {
	
	int i,j,k;

	FILE *file1;
	file1 = fopen("Part_Matrix.txt", "w");
	for (k=0; k<newcount; k++) {
		fprintf(file1,"%7.0f\t%9.6f\t%8.6f\t%8.6f\t%9.6f\t%9.6f\t%9.6f\n",Part_Matrix[k][0],
				Part_Matrix[k][1],Part_Matrix[k][2],Part_Matrix[k][3],
				Part_Matrix[k][4],Part_Matrix[k][5],Part_Matrix[k][6]);}
	fclose(file1);

	FILE *file2;
	file2 = fopen("PS_Phi_Final.txt", "w");
	for (j=0; j<=Ly; j++) {
		for (i=0; i<=Lx; i++) {
			fprintf(file2, "%9.6f\t",phi[i][j]);}
		fprintf(file2, "\n");}
	fclose(file2);
	
	FILE *file3;
	file3 = fopen("PS_ni_Final.txt", "w");
	for (j=0; j<=Ly; j++) {
		for (i=0; i<=Lx; i++) {
			fprintf(file3, "%8.6f\t",norm_inode[i][j]);}
		fprintf(file3, "\n");}
	fclose(file3);
	
	FILE *file4;
	file4 = fopen("Particle_Count.txt", "w");
	for (k=0; k<=finaltime; k++) {
		fprintf(file4,"%8d\n",pct[k]);}
	fclose(file4);
}

double nl_phi_solver(double ni_phi[]) {
	
	int i,j,k,m,n; int counter = 0;
	double sum_mp,sum_norm; double resid_norm = 1,tol = 1e-3;
	
	double Bmatrix[AMIN];
	double phi_new[AMIN] = {0},phi_old[AMIN] = {0};
	double alpha[AMIN],beta[AMIN],Matrix_Product[AMIN],resid[AMIN];
	
	while (resid_norm > tol) {
		//	Update Phi, Alpha, Beta
		for (k=0; k<nmax; k++) {
			phi_old[k] = phi_new[k];
			alpha[k] = 4 + delta2*exp(phi_old[k]);
			beta[k] = (1 - phi_old[k])*exp(phi_old[k]) - ni_phi[k];}
		//	Update Amatrix & Bmatrix
		for (j=0; j<jmax; j++) {
			for (i=0; i<imax; i++) {
				k = imax*j + i;
				Amatrix[k][k] = -alpha[k];
				if (i == 0) {Bmatrix[k] = delta2*beta[k] - LeftBC;}
				else if (i == imax-1) {Bmatrix[k] = delta2*beta[k] - RightBC;}
				else {Bmatrix[k] = delta2*beta[k];}}}
		//	Implement SOR Solver
		for (m=0; m<nmax; m++) {
			sum_mp = 0; sum_norm = 0;
			for (n=0; n<nmax; n++) {
				sum_mp += Amatrix[m][n]*phi_new[n];}
			Matrix_Product[m] = sum_mp;}
		for (k=0; k<nmax; k++) {
			resid[k] = Matrix_Product[k] - Bmatrix[k];
			phi_new[k] = phi_old[k] + omega*resid[k]/alpha[k];
			sum_norm += resid[k]*resid[k];}
		resid_norm = sqrt(sum_norm);
		counter++;}
	
	for (k=0; k<nmax; k++) {
		ni_phi[k] = phi_new[k];}
	
	printf("\nTolerance:%13.6f\n",tol);
	printf("Residual:%14.6f\n",resid_norm);
	printf("Counter:%15d\n",counter);
	
	return(1);
}

double rng(double num) {
	
	int MAXVALUE = 1e3;
	double nMAXVALUE;
	nMAXVALUE = (double)1/MAXVALUE;
	
	int value;
	value = rand() % (MAXVALUE + 1);
	if (value > MAXVALUE) {value = 0;}
	
	num = value*nMAXVALUE;
	
	return (num);
}
