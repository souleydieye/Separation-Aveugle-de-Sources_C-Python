///////////////////////////////////////
// separation_v0.c                   //
// Separation of sources with SOBI   //
// Test function: matrix A is known  //
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utilsFunctions.h"
#include "SOBI.h"

void separation_v0(double * input_1, double * input_2, double * output_1, double * output_2,
                        int nb_sample, int time_delay, int step){
  double normalisation;
  double A11, A12, A21, A22, D;
  int i;
  double **A_chapo = malloc(step * sizeof(*A_chapo));
  double **R_S1S1 = malloc(step * sizeof(*R_S1S1));
  double **R_S2S2 = malloc(step * sizeof(*R_S2S2));
  double **R_S1S2 = malloc(step * sizeof(*R_S1S2));

  normalisation = sqrt(2/(std(input_1,nb_sample)*std(input_1,nb_sample)+std(input_2,nb_sample)*std(input_2,nb_sample)));
  for (int i=0; i<nb_sample; i++){
    input_1[i] *= normalisation;
    input_2[i] *= normalisation;
  }
  normalizedR_S1S2(R_S1S1 ,input_1, input_1, nb_sample, step);
  normalizedR_S1S2(R_S2S2 ,input_2, input_2, nb_sample, step);
  normalizedR_S1S2(R_S1S2,input_1, input_2, nb_sample, step);

  transformationMatrix(A_chapo, R_S1S1, R_S2S2, R_S1S2, nb_sample, step);

  A11 =  A_chapo[0][0];
  A12 = A_chapo[0][1];
  A21 = A_chapo[1][0];
  A22 = A_chapo[1][1];

  D   = A11*A22-A12*A21;

  printf("\n A  = [[%f, %f]\n      [%f, %f]]\n\n", A11, A12, A21, A22);

  for (i = 0; i<nb_sample; i++) {
    output_1[i] =  (A22*input_1[i] - A12*input_2[i])/D;
    output_2[i] =  (-A21*input_1[i] + A11*input_2[i])/D;
  }

  for (int j = 0; j < step; j++) {
        free(R_S1S1[j]);
        free(R_S2S2[j]);
        free(R_S1S2[j]);
        free(A_chapo[j]);
  }
  free(R_S1S1);
  free(R_S2S2);
  free(R_S1S2);
  free(A_chapo);

}
