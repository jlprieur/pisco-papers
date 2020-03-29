/************************************************************************
* To compute O-C: delta_rho, delta_theta, delta_position
*
* Author: JLP
* Version: 26/04/2012 
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef PI
#define PI 3.14159
#endif

#define ABS(a) ((a < 0) ? (-(a)) : (a))

int main(int argc, char *argv[])
{
double rho_O, theta_O, rho_C, theta_C;
double drho, dtheta, dtheta_rad, dposi;
int nval1, nval2;
/* Check the number of arguments: 
*/
 fprintf(stderr, "argc=%d \n", argc);
if(argc > 1 && argc != 3) {
 fprintf(stderr, "Error: syntax is compu_resid rho_O,theta_O rho_C,theta_C \n");
 exit(-1);
 } else if (argc == 3) {
 nval1 = sscanf(argv[1], "%lf,%lf", &rho_O, &theta_O);
 nval2 = sscanf(argv[2], "%lf,%lf", &rho_C, &theta_C);
 } else {
 printf("rho_O (arcsec), theta_O (degr): ");
 nval1 = scanf(argv[1], "%lf,%lf", &rho_O, &theta_O);
 printf("rho_C (arcsec), theta_C (degr): ");
 nval2 = scanf(argv[2], "%lf,%lf", &rho_C, &theta_C);
 }

 if(nval1 != 2 || nval2 != 2) {
    fprintf(stderr, "Fatal error reading rho_O,theta_O rho_C,theta_C \n");
    fprintf(stderr, "argv1=>%s< argv2=>%s< \n", argv[1], argv[2]);
    exit(-1);
   }

printf("OK: rho_O=%f theta_O=%f rho_C=%f theta_C=%f\n",
       rho_O, theta_O, rho_C, theta_C);
drho = rho_O - rho_C;
dtheta = theta_O - theta_C;
dtheta_rad = dtheta * PI / 180.;
dposi = sqrt(drho * drho + rho_O*rho_O * dtheta_rad * dtheta_rad);

printf("Delta rho_(O-C)=%.3f (arcsec)   Delta_theta_(O-C)=%.2f (deg) \n",
       drho, dtheta);
printf("Delta_posi_(O-C)=%.3f (arcsec) Delta_posi/rho_O=%.2f (%)\n", dposi, 100. * dposi /(ABS(rho_O)+0.0001));

return(0);
}
