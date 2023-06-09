#include "bits/stdc++.h"
#include "deviates.h" // include this file last to avoid errors

using namespace std;

// parameters for EC-EC dynamics
// Timescale 5e4
const int N = 101;    // Number of grid points
const double L = 1.0; // Length of domain
const double dx = 1; // Grid spacing (dx = dy)
const double t_lambda = 1.0;
const double dt = 0.0001*t_lambda;    // Time step
const double D_phi = 5.0e-6;   // Transport Coefficient
const double kb = 5.0e-1; // birth rate
const double a = 1.0; // noise term (1 turn on & 0 turn off)

const int time_steps = 10; // timesteps

// parameters for colored noise
const double vs = 35.0; // self propulsion speed
const double Tp = 2.0e-1; // persistent time

double u[N][N];
double un[N][N];
// spatial noise values previous time
double xi[N][N];
double dxi[N][N];
double div_noise[N][N];
double brown[N][N];

// define random variables
Normaldev white_noise(0,sqrt(dt),time(0));
Normaldev stationary(0,vs/sqrt(2.0),time(0)); // equilibrium distribution of OU process
Ran ran(time(0));

// Periodic boundary conditions
void bc(){
	for(int i=0;i<N;i++){
		u[i][0] = u[i][N-1];
		u[0][i] = u[N-1][i];

		xi[i][0] = xi[i][N-1];
		xi[0][i] = xi[N-1][i];

		brown[i][0] = brown[i][N-1];
		brown[0][i] = brown[N-1][i];		
	}

	// make all corners equal
	u[N-1][0] = u[0][0];
	u[0][N-1] = u[0][0];
	u[N-1][N-1] = u[0][0];

	xi[N-1][0] = xi[0][0];
	xi[0][N-1] = xi[0][0];
	xi[N-1][N-1] = xi[0][0];

	brown[N-1][0] = brown[0][0];
	brown[0][N-1] = brown[0][0];
	brown[N-1][N-1] = brown[0][0];
}

// Initial condition
void init(){
	// initialize noise and set bc
	for(int x=0;x<N;x++){
		for(int y=0;y<N;y++){
			xi[x][y] = stationary.dev();
			brown[x][y] = white_noise.dev();
		}
	}

	bc();

	// initialize to uniform density
	for(int x=0;x<N;x++){
		for(int y=0;y<N;y++){
			//u[i] = exp(-pow((i*dx-L/2), 2) / 0.1); // Gaussian pulse centered in the middle of domain
			u[x][y] = 50.0/(1.0e4);//+ (1/100.0)*(ran.doub()-0.5); // uniform density (confluent ie set equal to 1/area)
			un[x][y] = u[x][y];

			div_noise[x][y] = (xi[x][y] + brown[x][y])*sqrt(un[x][y]);
		}
	}
}

void spde_solver(){
	// create files to output data
	ofstream pde("ec_ec_spde_2D.txt");

	// initialize
	init();
	bc();

	// output to file
	for(int x=0;x<N;x++){
		for(int y=0;y<N;y++){
			pde << u[x][y] << " ";
		}
		pde << endl;
	}

	// integrate with time
	for(int t=0;t<time_steps;t++){
		for(int y=0;y<N;y++){
			for(int x=0;x<N;x++){
				if((x==0) && (y==0)){ // first corner
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[y+1][x]+un[y][x+1]+un[y-1][x]+un[y][x-1]-4*(un[y][x])) + (dt/dx)*(div_noise[y][x+1]+div_noise[y+1][x]-2*div_noise[y][x]);
				} else if((x>0) && (x<N-1) && (y==0)){
					// First Edge
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[y+1][x]+un[y][x+1]+un[N-2][x]+un[y][x-1]-4*(un[y][x])) + (dt/dx)*(div_noise[y][x+1]+div_noise[y+1][x]-2*div_noise[y][x]);
				} else if((x==N-1) && (y==0)){
					// Second corner
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[y+1][x]+un[y][1]+un[N-2][x]+un[y][x-1]-4*(un[y][x])) + (dt/dx)*(div_noise[y][1]+div_noise[y+1][x]-2*div_noise[y][x]);
				} else if((x==0) && (y<N-1) && (y>0)){
					// Second Edge
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[y+1][x]+un[y][x+1]+un[y-1][x]+un[y][N-2]-4*(un[y][x])) + (dt/dx)*(div_noise[y][x+1]+div_noise[y+1][x]-2*div_noise[y][x]);
				} else if((x==0) && (y==N-1)){
					// Third corner
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[1][x]+un[y][x+1]+un[y-1][x]+un[y][N-2]-4*(un[y][x])) + (dt/dx)*(div_noise[y][x+1]+div_noise[1][x]-2*div_noise[y][x]);
				} else if((x==N-1) && (y<N-1) && (y>0)){
					// Third Edge
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[y+1][x]+un[y][1]+un[y-1][x]+un[y][x-1]-4*(un[y][x])) + (dt/dx)*(div_noise[y][1]+div_noise[y+1][x]-2*div_noise[y][x]);
				} else if((x==N-1) && (y==N-1)){
					// Fourth corner
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[1][x]+un[y][1]+un[y-1][x]+un[y][x-1]-4*(un[y][x])) + (dt/dx)*(div_noise[y][1]+div_noise[1][x]-2*div_noise[y][x]);
				} else if((x>0) && (y==N-1) && (x<N-1)){
					// Fourth edge
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[1][x]+un[y][x+1]+un[y-1][x]+un[y][x-1]-4*(un[y][x])) + (dt/dx)*(div_noise[y][x+1]+div_noise[1][x]-2*div_noise[y][x]);
				} else{
					u[y][x] = un[y][x] + ((D_phi*dt)/pow(dx,2))*(un[y+1][x]+un[y][x+1]+un[y-1][x]+un[y][x-1]-4*(un[y][x])) + (dt/dx)*(div_noise[y][x+1]+div_noise[y+1][x]-2*div_noise[y][x]);
				}
			}
		}

		bc();

		// calculate noise terms for next time steps
		for(int x=0;x<N;x++){
			for(int y=0;y<N;y++){
				un[x][y] = u[x][y];

				dxi[x][y] = -(xi[x][y]/Tp)*dt + (vs/sqrt(Tp))*white_noise.dev();
				xi[x][y] += dxi[x][y];
				brown[x][y] = white_noise.dev();
			}
		}

		bc();

		for(int x=0;x<N;x++){
			for(int y=0;y<N;y++){
				div_noise[x][y] = (xi[x][y] + brown[x][y])*sqrt(un[x][y]);
			}
		}

		// output to file
		for(int x=0;x<N;x++){
			for(int y=0;y<N;y++){
				pde << u[x][y] << " ";
			}
			pde << endl;
		}
	}
}


int main(){
	spde_solver();

	return 0;
}