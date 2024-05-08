#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string> 
#include <fstream>
#include <time.h>
//#define SNR_dB 10.0      // Signal-to-noise ratio in dB
//#define NUM_BITS 5000000
//#define EbN0 (pow(10, (SNR_dB / 10.0)))
#define PI 3.14159265358979323846
#define Ts 0.00003
#define Wm 2*PI*200 
#define N 32.0
#define N0 8.0

using namespace std;

int main() {
	
	int SNR_dB;      // Signal-to-noise ratio in dB
	int NUM_BITS = 5000000;
    double EbN0;
    int j, Tx, Rx, Sn, Sn_1 = 1;
    int num_errors = 0;
    double hi_t_real, hi_t_complex, hq_t_real, hq_t_complex;
    double rn, rn_1, hi_t, yn, yn_1, m, sim_ber, theo_ber ;
    unsigned seed;
    seed = (unsigned)time(NULL); 
    srand(seed); 
	ofstream file("dbpsk_data.csv");
	//file << "0" << "," << "2.82843" << "," << "0" <<","<< "2.82843" << "," << "0" <<endl;
	rn_1 = sqrt(2.82843*2.82843);
	yn_1 = rn_1*Sn_1 + gaussian_noise(EbN0);
	
	for(SNR_dB = 0; SNR_dB<=50; SNR_dB+=5){
	
		EbN0 = (double)(pow(10, (SNR_dB / 10.0)));
		
		for(j = 1; j<NUM_BITS; j++){
			
			Tx = rand() % 2;
			if(Tx==0) Sn = Sn_1 * -1 ;
			else Sn = Sn_1;
					
			hi_t_real = h_t_real(j*Ts);
			hq_t_complex = h_t_complex(j*Ts);
			rn = sqrt(hi_t_real*hi_t_real + hq_t_complex*hq_t_complex);
			
			yn = (double)rn*Sn + gaussian_noise(EbN0);
			m = yn*yn_1;
			
			if(m>=0) Rx = 1;
			else Rx = 0;
			
			if(Rx!=Tx) num_errors++;
			
			yn_1 = yn;
	    	rn_1 = rn;
	    	Sn_1 = Sn;	
		}
		
		sim_ber = (double)num_errors / (NUM_BITS);
		theo_ber = 0.5*(1- j0(Wm*Ts)/(1+1/EbN0) );

		printf("Sim_Ber: %f\n", log10(sim_ber));
		printf("Theo_Ber: %f\n", log10(theo_ber));
		printf("Errors: %d\n", num_errors);
		//write to csv
    	file << SNR_dB << "," << sim_ber <<endl;
    	yn_1 = yn;
    	rn_1 = rn;
    	Sn_1 = Sn;
    	num_errors = 0;
	}
	file.close();
    return 0;
}

// Function to generate Gaussian-distributed random noise
double gaussian_noise(double EbN0) {
    double u1, u2, n;
    u1 = (double)(rand()) / (RAND_MAX + 1.0) ;
    u2 = (double)(rand()) / (RAND_MAX + 1.0) ;    
    n = (double) sqrt(1.0/(2.0*EbN0) ) * sqrt(-2.0 * log(u1)) * sin(2.0 * PI * u2);
    return n;
}

// compute alpha
double alpha(int i){
 	
 	double al = (2*PI*i)/N + PI/(2*N);
 	return al;
 }
 
double phi_nk(int i){
	
	double phi_n = (i*PI)/N0;
	return phi_n;
	
}
// compute real part
double h_t_real(double t){
	int i;
	double al,ht,sum = 0,cos_al;	
	for(i = 0; i<N0; i++){
		al = alpha(i);
		cos_al = cos(al);
		sum += cos(Wm*t*cos_al);
	}
	ht = (double) (1.0/sqrt(N0))*sum;
	return ht;
}
// compute real part
double h_t_complex(double t){
	int i;
	double al,ht,sum = 0,sin_al;
	for(i = 0; i<N0; i++){
		al = alpha(i);
		sin_al = sin(al);
		sum += sin(Wm*t*sin_al);
	}
	ht = (double) (1.0/sqrt(N0))*sum;
	return ht;
}


