#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string> 
#include <fstream>
#define SNR_dB 10.0      // Signal-to-noise ratio in dB
#define EbN0 (pow(10, (SNR_dB / 10.0)))
#define PI 3.14159265358979323846
#define Ts 0.00003
#define Wm 2*PI*80
#define N 32.0
#define N0 8.0

using namespace std;


int main() {
    int j;
    int NUM_BITS = 5000000;
    int num_errors = 0;
    double hi_t_real, hi_t_complex, hq_t_real, hq_t_complex;
    double real_part, complex_part, rn,hi_t;

	ofstream file("dbpsk_data.csv");
	
	for(j = 0; j<NUM_BITS; j++){
		hi_t_real = h_t_real(j*Ts);
		hq_t_complex = h_t_complex(j*Ts);

		rn = sqrt(hi_t_real*hi_t_real + hq_t_complex*hq_t_complex);
		
		double phi = atan2(hq_t_complex,hi_t_real);
		//write to csv
    	//file << real_part << "," << complex_part << "\n"<< endl;
    	file << j*Ts << "," << rn<< "," << phi <<","<< hi_t_real << "," << hq_t_complex <<endl;	
	}
	file.close();
    return 0;
}

// Function to generate Gaussian-distributed random noise
double gaussian_noise() {
    
    double u1 = rand() / (RAND_MAX + 1.0);
    double u2 = rand() / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
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
	ht = (double) (1/sqrt(N0))*sum;
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
	ht = (double) (1/sqrt(N0))*sum;
	return ht;
}
