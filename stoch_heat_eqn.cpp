#include "bits/stdc++.h"
#include "deviates.h" // include this file last to avoid errors

using namespace std;

// parameters for EC-EC dynamics
// Timescale 5e4
const int N = 101;    // Number of grid points
const double L = 1.0; // Length of domain
const double dx = 1; // Grid spacing
const double t_lambda = 1.0;
const double dt = 0.001*t_lambda;    // Time step
const double D_phi = 5.0e-6;   // Transport Coefficient
const double kb = 0.5; // birth rate
const double a = 1.0; // noise term (1 turn on & 0 turn off)

const int iter = 10000; // timesteps

// parameters for colored noise
const double vs = 35.0; // self propulsion speed
const double Tp = 0.20; // persistent time


double u[N]; // current values
double un[N]; // value at previous time step
// spatial noise values previous time
double xi[N];
double dxi[N];
double div_noise[N];

// define noise variables
Normaldev white_noise(0,sqrt(dt),time(0));
Normaldev stationary(0,vs/sqrt(2.0),time(0));
Ran ran(time(0));

// Initial condition
void init(){
	// initialize to uniform density
	for(int i=0;i<N;i++){
		xi[i] = stationary.dev();

		//u[i] = exp(-pow((i*dx-L/2), 2) / 0.1); // Gaussian pulse centered in the middle of domain
		u[i] = 10.0/100.0 + (1/100.0)*(ran.doub()-0.5); // uniform density
		un[i] = u[i];

		div_noise[i] = (xi[i] + white_noise.dev())*sqrt(un[i]);

	}
}

// Periodic boundary conditions
void bc(){
	u[0] = u[N-1];
}

// colored noise
// autocorrelation follows Ornstein-Uhlenbeck process
// dx = -kx dt + sqrt(D) dW

/*void colored_noise(){
	// Create timeseries starting from sampled stationary distribution to get xi(t)
	Normaldev white_noise(0,sqrt(dt),time(0));
	Normaldev stationary(0,vs/sqrt(2.0),time(0));
	double dxi, xi;
	ofstream colored("colored_noise.txt");

	colored << 0 << endl;

	xi = stationary.dev();

	//colored << xi << endl;

	for(int i=0;i<10;i++){
		dxi = -(xi/Tp)*dt + (vs/sqrt(2.0))*white_noise.dev(); // Langevin equation
		cout << dxi << endl;
		xi = xi + dxi;
		//colored << xi << endl;
	}
}*/

// simulation to solve stochastic diffusion equation
void stoch_heat_eqn_solver(){
	ofstream fout("pde_data.txt");
	ofstream fout2("x_axis.txt");
	ofstream fout3("noise_term.txt");
	ofstream ou("ou_process.txt");

	ou << 0 << endl;

	// x axis data
	fout2 << "x axis" << endl;
	for(int i=0;i<N;i++){
		fout2 << i*dx << endl;
		fout << i << " ";
		fout3 << i << " ";
	}
	fout3 << endl;
	fout << endl;

	init();
	bc();

	//ou << xi[0] << endl;

	// store initial data
	for(int i=0;i<N;i++){
		fout << u[i] << " ";
	}
	fout << endl;
	
	//Normaldev rand(0.0,sqrt(2.0*alpha*dt),time(0));

	// integrate in time 10000 iterations
	for(int j=0;j<iter;j++){
		// finite difference scheme

		// include noise term in this loop

		for(int i=0;i<N;i++){
			if(i==0){
				u[i] = un[i] + D_phi*(dt/pow(dx,2))*(un[i+1] - 2*un[i] + un[N-1]) + (kb*un[i]*dt) + a*(div_noise[i+1] - div_noise[i])*(dt/dx);
			} else if(i==N-1){
				u[i] = un[i] + D_phi*(dt/pow(dx,2))*(un[0] - 2*un[i] + un[i-1]) + (kb*un[i]*dt) + a*(div_noise[0] - div_noise[i])*(dt/dx);
			}else{
				u[i] = un[i] + D_phi*(dt/pow(dx,2))*(un[i+1] - 2*un[i] + un[i-1]) + (kb*un[i]*dt) + a*(div_noise[i+1] - div_noise[i])*(dt/dx);
			}

		}

		bc();

		// output to file
		for(int i=0;i<N;i++){
			fout << u[i] << " ";
			fout3 << div_noise[i] << " ";
		}
		fout3 << endl;
		fout << endl;

		for(int i=0;i<N;i++){
			un[i] = u[i];
			dxi[i] = -(xi[i]/Tp)*dt + (vs/sqrt(Tp))*white_noise.dev(); // for next time step
			xi[i] += dxi[i];
			div_noise[i] = (xi[i] + white_noise.dev())*sqrt(un[i]);
		}
		if (j==1000){
			for (int i=0;i<N;i++){
				ou << xi[i] << endl;
			}
		}
		//ou << xi[0] << endl; // store ou process values for first coordinate
	}
}

int main(){
	//colored_noise();
	stoch_heat_eqn_solver();
	/*Normaldev rand(0,sqrt(2*0.0001),time(0));
	ofstream fout3("gaussian_noise.txt");
	double noise[100];
	fout3 << 0 << endl;
	for(int i=0;i<100;i++){
		noise[i] = rand.dev();
		fout3 << noise[i] << endl;
	}*/

	return 0;
}