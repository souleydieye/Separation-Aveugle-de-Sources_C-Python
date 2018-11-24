// -----------------------------------------------------------------//
// Code property:                                                   //
// ==============                                                   //
// TELECOM BRETAGNE						    //
// Dpt. Signal et Communications				    //
// Technopole Brest-Iroise					    //
// CS 83818 - 29238 Brest Cedex 3 - France			    //
//                                                                  //
// General info:                                                    //
// =============                                                    //
// program: math_functions.c                                        //
// last upate: 14/09/2016                                           //
// info: thierry.chonavel@telecom-bretagne.eu                       //
//       thierry.legall1@telecom-bretagne.eu                        //
//                                                                  //
// Code info:  	                                                    //
// ==========                                                       //
// A bunch of math functions for program Sobi.c                     //
//------------------------------------------------------------------//
#include <math.h>

// empirical mean of vector entries
double mean(double * input, int nb_samples) {
	double sum = 0.0;
	for (int i=0; i<nb_samples; i++) sum += input[i];
  return sum/nb_samples;
}
// center vector by substracting mean
void center(double * input, int nb_samples) {
	double m = mean(input, nb_samples);
	for (int i=0; i<nb_samples; i++) input[i] -= m;
}
// dot product of two vectors
double dot_product(double * input_1, double * input_2, int nb_samples) {
	double dot_prod = 0.0;
	for (int i=0; i<nb_samples; i++) dot_prod += input_1[i]*input_2[i];
	return dot_prod;
}
// correlation of two vectors
double correlation(double * input_1, double * input_2, int nb_samples) {
	double s = 0;
	double m1 = mean(input_1, nb_samples);
	double m2 = mean(input_2, nb_samples);
	for (int i=0; i<nb_samples; i++) s+= (input_1[i]-m1)*(input_2[i]-m2);
	return s/nb_samples;
}
// squared vector norm
double squared_norm(double * input, int nb_samples) {
	return dot_product(input, input, nb_samples);
}
// vector norm
double norm(double * input, int nb_samples) {
	return sqrt(squared_norm(input, nb_samples));
}
// vector standard deviation
double std(double * input, int nb_samples) {
	double s = 0;
	double m = mean(input,nb_samples);
	for (int i=0; i<nb_samples; i++) s+= (input[i]-m)*(input[i]-m);
	return sqrt(s/nb_samples);
}


void normalizedR_S1S2(double **matrixRs1s2, double * input_1, double * input_2, int nb_samples, int step)
{
    double walk;
    walk = nb_samples / step;


    for (int i = 0; i < step; i++){
        matrixRs1s2[i] = malloc(step * sizeof(**matrixRs1s2));
	}

    for (int i = 0; i < step; i++){
        for (int j=0; j< step; j++){
            matrixRs1s2[i][j] = 0.0;
        }
    }

    for (int i = 0; i < walk; i++){
		for (int j = 0; j < step; j++){
			for (int k =0; k < step; k++){
                matrixRs1s2[j][k] += input_1[j + i*step] * input_2[k + i*step];
            }
        }
    }

    for (int i = 0; i < step; i++){
			for (int j =0; j< step; j++){
                matrixRs1s2[i][j] /= walk ;
            }
        }

}

double normalizedTrace(double **matrixRs1s2,int step)
{
    double trace=0.0;


    for (int i = 0; i < step; i++){
        trace += matrixRs1s2[i][i];
    }
    return trace/step;
}

double sum(double **matrixRs1s2,int step)
{
    double sum=0.0;

    for (int i = 0; i < step; i++) {
        for (int j = 0; j < step; j++) {
            sum += matrixRs1s2[i][j];
        }
    }

    return sum;
}

double off(double **matrixRs1s2,int step)
{
	double off=0.0;

    off = (sum(matrixRs1s2,step)-step*normalizedTrace(matrixRs1s2,step)) / ( step * (step-1) );
    return off;
}
