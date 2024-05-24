/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Program calib_rho 
* To compute the scale for speckle interferometry
* using autocorrelations images of Deneb with a multi-slit mask
*
* JLP 
* Version 17-10-2011
-------------------------------------------------------------------*/
#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
// #include <jlp_ftoc.h>
#include "jlp_numeric.h"
#include "jlp_fitsio.h"

#define N_FILTERS 5

/*
* PISCO filters with ICCD and Deneb (A2 I type-star, #120):
*
* Filter B    Cent. wavel: 4521 Bandwidth: 439
* Filter OIII Cent. wavel: 5010 Bandwidth: 110
* Filter V    Cent. wavel: 5295 Bandwidth: 558
* Filter R    Cent. wavel: 6406 Bandwidth: 697
* Filter RL   Cent. wavel: 7387 Bandwidth: 687
* Filter IR   Cent. wavel: 8420 Bandwidth: 643 (not used for calibration)
*/

int main(int argc, char *argv[])
{
char in_name[64], filter_name[20];
double aa_mm, error_aa_mm, lambda[N_FILTERS], dlambda[N_FILTERS]; 
double scale, sum, sumsq, mean_scale, error_scale, temperature;
double rho, error_rho, theta_arcsec, delta_aa;
char filt_name[N_FILTERS][20];
FILE *fp_in;
register int i;

aa_mm = 87.96;
error_aa_mm = 0.02;

lambda[0] = 452.1;
dlambda[0] = 43.9;
strcpy(filt_name[0], "B");

lambda[1] = 501.0;
dlambda[1] = 11.0;
strcpy(filt_name[1], "OIII");

lambda[2] = 529.5;
dlambda[2] = 55.8;
strcpy(filt_name[2], "V");

lambda[3] = 640.6;
dlambda[3] = 69.7;
strcpy(filt_name[3], "R");

lambda[4] = 738.7;
dlambda[4] = 68.7;
strcpy(filt_name[4], "RL");

/*
lambda[5] = 842.0;
dlambda[5] = 64.3;
strcpy(filt_name[5], "IR");
*/

printf(" Program calib_rho to compute the scale \n");
printf(" JLP Version 12-05-2023 \n");

/* One parameters only is allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 3 && argc != 2 && argc != 1)
  {
  printf(" Syntax: calib_rho in_file temperature_deg_C \n"); 
  printf(" in_file: file with autocorrelation measurements\n");
  printf(" Fatal: Syntax error: argc=%d\n",argc);
  exit(-1);
  }

/* Input of parameters with the command line: */
if (argc == 3 )
 { 
  strcpy(in_name, argv[1]);
  sscanf(argv[2], "%lf", &temperature);
 }
else if (argc == 2 )
 { 
  strcpy(in_name, argv[1]);
  temperature = 20.;
 }
/* Interactive input of parameters: */
else
 { 
  printf(" Input file := "); scanf("%s", in_name);
  printf(" Temperature (deg. C):= "); scanf("%lf", &temperature);
 }

if((fp_in = fopen(in_name, "r")) == NULL) {
 fprintf(stderr, "Fatal error opening input file: >%s<\n", in_name);
 return(-1);
 }

delta_aa = aa_mm * (temperature - 20) * 24.e-6;
aa_mm += delta_aa;
printf(" T=%.2f aa_mm=%.2f Delta=%.2f\n", temperature, aa_mm, delta_aa);
sum = 0;
sumsq = 0.;
for(i = 0; i < N_FILTERS; i++) {
  fscanf(fp_in,"%lf %lf %s\n", &rho, &error_rho, filter_name);
  theta_arcsec = (lambda[i] * 1.e-9 / (aa_mm * 1.e-3)) * 180. * 3600. / PI;
  scale = theta_arcsec / rho;
  printf("theta = %f arcsec rho = %f pixels\n", theta_arcsec, rho);
  error_scale = SQUARE(error_rho / rho) + SQUARE(dlambda[i] / (2. * lambda[i]))
                + SQUARE(error_aa_mm / aa_mm);
  error_scale = sqrt(error_scale);
  sum += scale;
  sumsq += SQUARE(scale);
  printf("scale for %s: %f +/- %f (= %s ?)\n", filt_name[i], scale, 
          error_scale, filter_name);
  }

mean_scale = sum / (double)N_FILTERS;
error_scale = sqrt(sumsq / (double)N_FILTERS - SQUARE(mean_scale));
printf("Mean scale = %f +/- %f\n", mean_scale, error_scale);

fclose(fp_in);

return(0);
}
