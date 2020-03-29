/************************************************************************
* To compute magnitude of components
* from Mag_total and Delta_mag
*
* Author: JLP
* Version: 04/05/2012 
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ABS(a) ((a < 0) ? (-(a)) : (a))

int main(int argc, char *argv[])
{
double total_mag, dmag, mt, mag1, mag2; 
int nval, iter;
/* Check the number of arguments: 
*/
 fprintf(stderr, "argc=%d \n", argc);
if(argc > 1 && argc != 2) {
 fprintf(stderr, "Error: syntax is compu_mag Total_mag,Delta_mag \n");
 exit(-1);
 } else if (argc == 2) {
 nval = sscanf(argv[1], "%lf,%lf", &total_mag, &dmag);
 } else {
 printf("total_mag, delta_mag: ");
 nval = scanf(argv[1], "%lf,%lf", &total_mag, &dmag);
 }

 if(nval != 2) {
    fprintf(stderr, "Fatal error reading Total_mag, Delta_mag \n");
    fprintf(stderr, "argv1=>%s< \n", argv[1]);
    exit(-1);
   }

printf("OK: total_mag=%f dmag=%f\n", total_mag, dmag);

/*
mag = -2.5 log10 (Flux)
Flux = 10**(-0.4 * mag)
total_mag = -2.5 log10 [10**(-0.4 * mag1) + 10**(-0.4 * mag2)]
*/

/* Iterative resolution: */
mag1 = total_mag;
for(iter = 0; iter < 100000; iter++) {
mag2 = mag1 + dmag;
mt = -2.5 * log10(pow(10.,(-0.4 * mag1)) + pow(10.,(-0.4 * mag2)));
/*
printf("mt = %.4f total_mag = %.4f\n", mt, total_mag); 
printf("mag1 = %.4f mag2 = %.4f\n", mag1, mag2); 
*/
if(ABS(mt - total_mag) < 0.001) break;
mag1 -= (mt - total_mag)/2.;
}

printf("iter= %d, mt = %.4f, total_mag = %.4f\n", iter, mt, total_mag); 
printf("mag_1=%.3f mag_2=%.3f \n", mag1, mag2);

return(0);
}
