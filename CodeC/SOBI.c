#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utilsFunctions.h"


void transformationMatrix(double **A, double **normalizedR_S1S1, double **normalizedR_S2S2, double **normalizedR_S1S2, int nb_sample, int step)
{
	double F1,F2,F12,T1,T2,T12,alpha,beta,gamma,d1,d2;

	for (int i = 0; i < step; i++){
        A[i] = malloc(step * sizeof(**A));
	}

    for (int i = 0; i < step; i++){
        for (int j=0; j<step; j++){
            A[i][j] = 0.0;
        }
    }

    T1 = normalizedTrace(normalizedR_S1S1, step);
    T2 = normalizedTrace(normalizedR_S2S2, step);
    T12 = normalizedTrace(normalizedR_S1S2, step);

    F1 = off(normalizedR_S1S1, step);
    F2 = off(normalizedR_S2S2, step);
    F12 = off(normalizedR_S1S2, step);

    alpha = 2.0*F12*T12 -(F1*T2 + F2*T1);
    beta = 2.0*(pow(T12,2) -T1*T2);
    gamma = sqrt( pow((F1*T2 -F2*T1),2) + 4.0*(F12*T2 - T12*F2)*(F12*T1 - T12*F1));

    d1 = alpha - gamma;
    d2 = alpha + gamma;


    A[0][0] = beta*F1 - T1*d1;
    A[0][1] = beta*F12 -T12*d2;
    A[1][0] = beta*F12 -T12*d1;
    A[1][1] = beta*F2-T2*d2;

}


























