// -----------------------------------------------------------------//
// header for math_functions.c                                      //			                                      //
//------------------------------------------------------------------//
#include <math.h>

// empirical mean of vector entries
double mean(double * input, int nb_samples);
// center vector by substracting mean
void center(double * input, int nb_samples);
// dot product of two vectors
double dot_product(double * input_1, double * input_2, int nb_samples);
// correlation of two vectors
double correlation(double * input_1, double * input_2, int nb_samples);
// vector norm
double norm(double * input, int nb_samples);
// vector standard deviation
double std(double * input, int nb_samples);
//Correlation matrix between input_1 and input_2
void normalizedR_S1S2(double **matrixRs1s2, double * input_1, double * input_2, int nb_sample, int step);
//Normalized trace
double normalizedTrace(double **matrixRs1s2,int step);
// return the sum of all coefficients of our matrix
double sum(double **matrixRs1s2,int step);
//the off by definition
double off(double **matrixRs1s2,int step);

