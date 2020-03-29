/************************************************************************
* "residuals_utils.c"
*
* To compute the residuals of measurements of binary stars
* from known orbits.
*
* From o-c.for (version of 2008)
* written by Marco SCARDIA - Osservatorio Astronomico di Brera-Merate
*
*
* JLP 
* Version 08/04/2009
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "residuals_utils.h"

#define DEBUG
/*
#define DEBUG_1
*/

/*************************************************************
* Compute the ephemerids corresponding to the observation epoch
* of a multiple system
*
* OUTPUT:
*  theta_c = computed position angle (degrees)
*  rho_c= computed separation angle (arcseconds) 
*
*************************************************************/
int compute_ephemerid_of_multiple_system(int nber_of_orbits,
                      double *Omega_node, double *omega_peri, double *i_incl, 
                      double *e_eccent, double *T_periastron, double *Period, 
                      double *a_smaxis, double *mean_motion, double epoch_o, 
                      double *c_tolerance, double *theta_c, double *rho_c)
{
double wtheta, wrho;
double wx, wy;
register int k;

/* Compute the ephemerids corresponding to the observation epoch: */
wx = 0.;
wy = 0.;
for(k = 0; k < nber_of_orbits; k++) {
  compute_ephemerid(Omega_node[k], omega_peri[k], i_incl[k], e_eccent[k],
                    T_periastron[k], Period[k], a_smaxis[k], mean_motion[k],
                    epoch_o, c_tolerance[k], &wtheta, &wrho);
  wx += wrho * cos(wtheta * PI /180.);
  wy += wrho * sin(wtheta * PI /180.);
  }
*theta_c = atan2(wy,wx);
*rho_c = sqrt(wx*wx + wy*wy);

/* Conversion to degrees: */
 *theta_c *= (180./PI); 
/* Put the result into the interval [0,360] */
 if(*theta_c < 0.0) *theta_c += 360.0;
 if(*theta_c > 360.0) *theta_c -= 360.0;

return(0);
}
/*************************************************************
* Compute the ephemerids corresponding to the observation epoch: 
*
* INPUT:
*  Omega_node (radians), omega_peri (radians), i_incl (radians), 
*  e_eccent, T_periastron (years), Period (years),
*  a_smaxis (arcseconds) =  orbital elements
*  epoch_o = Epoch of the ephemerid
*  c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*
* OUTPUT:
*  theta_c = ANGOLO DI POSIZIONE CALCOLATO (degrees)
*  rho_c= DISTANZA ANGOLARE CALCOLATA
*
**************************************************************/
int compute_ephemerid(double Omega_node, double omega_peri, double i_incl, 
                      double e_eccent, double T_periastron, double Period, 
                      double a_smaxis, double mean_motion, double epoch_o, 
                      double c_tolerance, double *theta_c, double *rho_c)
{
double daa, aa, ab, pp, theta;
double ee, eps, cc, mean_anomaly, eccentric_anomaly, true_anomaly;
int isafe;
/* mean_anomaly = ANOMALIA MEDIA
* true_anomaly = ANOMALIA VERA
* eccentric_anomaly = ANOMALIA ECCENTRICA
*/
   daa = epoch_o - T_periastron;
/*   mean_anomaly = mean_motion * AMOD(daa, Period);
* real function that returns the value of daa modulo Period
* i.e.:
*/
   daa = daa - Period * (double)((int)(daa / Period));
   mean_anomaly = mean_motion * daa;
/* Initialize ee: */
   ee = mean_anomaly + e_eccent*sin(mean_anomaly)
             / sqrt(1.0 + e_eccent * (e_eccent - 2.0 * cos(mean_anomaly)));
/* Iterative resolution of Kepler's equation: */
   isafe = 0;
   do {
    eps = (mean_anomaly - ee) + e_eccent * sin(ee);
    cc = eps / (1.0 - e_eccent * cos(ee));
/* Correction of the current value of ee: */
    ee += cc;
    isafe++;
   } while((ABS(cc) > c_tolerance) && (isafe < 1000000));

   eccentric_anomaly = ee * 0.5;
   aa = sqrt((1.0 + e_eccent) / (1.0 - e_eccent)) * sin(eccentric_anomaly);
/* In FORTRAN: 
* ATAN2(y,x) returns the argument alpha of the complex number x + i y,
* expressed in radians, in the range [-PI, +PI],
* such that: tan(alpha) = y / x
* if x > 0  alpha = arctan(y/x)
* if x < 0  alpha = +PI + arctan(y/x) if y > 0
*       or  alpha = -PI + arctan(y/x) if y < 0
* if x = 0  alpha = PI/2 if y > 0 or alpha = -PI/2 if y < 0
*
* Same syntax in C: atan2(y,x) = arctan(y/x) = ATAN2(y,x) 
*/
   true_anomaly = 2.0 * atan2(aa, cos(eccentric_anomaly));
   pp = true_anomaly + omega_peri;
   if(pp > 2.*PI) pp -= 2.*PI;
   ab = cos(i_incl) * sin(pp);
   theta = atan2(ab, cos(pp));
   if(theta < 0.0) theta += 2.*PI;

/* Computed position angle: */
   *theta_c = theta + Omega_node;

/* Conversion to degrees: */
 *theta_c *= (180./PI); 
/* Put the result into the interval [0,360] */
 if(*theta_c < 0.0) *theta_c += 360.0;
 if(*theta_c > 360.0) *theta_c -= 360.0;


/* Computed separation angle: */
   *rho_c = a_smaxis * (1.0 - e_eccent * cos(ee)) * cos(pp) / cos(theta);
/* DEBUG:
printf("theta=%f theta_c=%f\n", theta, *theta_c);
printf("pp=%f ee=%f e_ccent=%f a_smaxis=%f rho_c=%f\n", pp, ee, e_eccent,
       a_smaxis, *rho_c);
*/

return(0);
}
/************************************************************
* Correction for precession
* using the formula of Armellini (1931) (or Couteau, 1978)
*
* INPUT : 
* theta: before correction (in radians)
* alpha, delta: coordinates of object (in radians)
* epoch_o: epoch of observation
* orbit_equinox: equinox used as a reference for computing the orbit 
*
* OUTPUT : 
* dtheta_precess: correction for precession (in radians)
*************************************************************/
int precession_correction(double *dtheta_precess, double alpha, double delta, 
                          double epoch_o, double orbit_equinox)
{

/* Precession correction of theta in arcseconds 
*/
  *dtheta_precess = -20.0 * (epoch_o - orbit_equinox) * sin(alpha) / cos(delta);

#ifdef DEBUG_1
  printf(" epoch=%f equinox=%f sin(alpha)=%f cos(delta)=%f\n", 
         epoch_o, orbit_equinox, sin(alpha), cos(delta));
  printf(" correction for precession: %f (arcsec) or %f (degrees)\n", 
         *dtheta_precess,  *dtheta_precess/3600.);
#endif

/* Conversion to radians: */
  *dtheta_precess *= DEGTORAD/3600.;

return(0);
}
/***************************************************************************
* read_orbital_elements
* Read input orbit file and retrieve the orbital parameters 
* in Marco's format if orbit_format = 1 
* or in OC6 format if orbit_format = 2
*
* INPUT:
* fp_in: pointer to the input file containing the measurements and the
*        orbital elements
* orbit_format: format of the input file (1=Marco, 2=OC6)
* nber_of_orbits: number of orbits to be found in input file
* 
* OUTPUT:
* arrays filled with the orbital parameters 
*
***************************************************************************/
int read_orbital_elements_from_file(char *orbit_infile, int orbit_format, 
              double *Omega_node, double *omega_peri, double *i_incl, 
              double *e_eccent, double *T_periastron, double *Period, 
              double *a_smaxis, double *mean_motion, double *orbit_equinox,
              int nber_of_orbits)
{
FILE *fp_in;
int nval, k, status;
char in_line[80];

if((fp_in = fopen(orbit_infile, "r")) == NULL) {
  fprintf(stderr, "Fatal error opening input file >%s<\n", orbit_infile);
  exit(-1);
  }

for(k = 0; k < nber_of_orbits; k++) {

/* Read next line: */
  while(!feof(fp_in)) {
    if(fgets(in_line, 80, fp_in)) {
      if(in_line[0] != '%') break;  
    }
  }

/* Read orbital elements from input file: */
/* Marco's format, i.e., orbit_format = -1 or orbit_format = 1 */
   if(orbit_format != -1 && orbit_format != 1) {
      fprintf(stderr, "read_orbital_elements_from_file/Fatal error: bad option: orbit_format = %d\n", 
              orbit_format);
      exit(-1);
      }
   nval = sscanf(in_line, "%lf %lf %lf %lf %lf %lf %lf %lf", 
                 &Omega_node[k], &omega_peri[k], &i_incl[k], &e_eccent[k], 
                 &T_periastron[k], &Period[k], &a_smaxis[k], &orbit_equinox[k]);
   if(nval != 8) {
     fprintf(stderr, "read_orbital_elements_from_file/Fatal error reading input file (nval=%d) \n %s \n", 
              nval, in_line);
     fclose(fp_in);
     return(-2);
    }

#ifdef DEBUG
 printf("nval=%d Omega_node=%.3f omega_peri=%.3f incl=%.3f e=%.3f T=%.3f P=%.3f a=%.3f Equinox=%.3f\n", 
        nval, Omega_node[k], omega_peri[k], i_incl[k], e_eccent[k], 
        T_periastron[k], Period[k], a_smaxis[k], orbit_equinox[k]);
#endif

/* Conversion to radians (necessary for ephemrids): */
  Omega_node[k] *= DEGTORAD;
  omega_peri[k] *= DEGTORAD;
  i_incl[k] *= DEGTORAD;
  mean_motion[k] = (360.0 / Period[k]) * DEGTORAD;
} /* EOF loop on k */

 fclose(fp_in);

if(k == nber_of_orbits) 
 status = 0;
else
 status = -1;

return(status);
}
/***************************************************************************
* read_orbital_elements_and_errors
* Read input orbit file and retrieve the orbital parameters 
* in Marco's format if orbit_format = 1 
* or in OC6 format if orbit_format = 2
*
* INPUT:
* fp_in: pointer to the input file containing the measurements and the
*        orbital elements
* orbit_format: format of the input file (1=Marco, 2=OC6)
* nber_of_orbits: number of orbits to be found in input file
* 
* OUTPUT:
* arrays filled with the orbital parameters 
*
***************************************************************************/
int read_orbital_elements_and_errors_from_file(char *orbit_infile, 
              int orbit_format, 
              double *Omega_node, double *omega_peri, double *i_incl, 
              double *e_eccent, double *T_periastron, double *Period, 
              double *a_smaxis, double *mean_motion, double *orbit_equinox,
              double *err_Omega_node, double *err_omega_peri, 
              double *err_i_incl, double *err_e_eccent, 
              double *err_T_periastron, double *err_Period, 
              double *err_a_smaxis, double *err_mean_motion, int nber_of_orbits)
{
FILE *fp_in;
int nval, i, k, status, expected_nval;
char in_line[80];

// Initialize errors to zero:
for(k = 0; k < nber_of_orbits; k++) {
  err_Omega_node[k] = 0.;
  err_omega_peri[k] = 0.;
  err_i_incl[k] = 0.; 
  err_e_eccent[k] = 0.;
  err_T_periastron[k] = 0.;
  err_Period[k] = 0.;
  err_a_smaxis[k] = 0.;
  err_mean_motion[k] = 0.;
  }

if((fp_in = fopen(orbit_infile, "r")) == NULL) {
  fprintf(stderr, "Fatal error opening input file >%s<\n", orbit_infile);
  exit(-1);
  }

k = -1;
for(i = 0; i < 2*nber_of_orbits; i++) {
/* Read next line: */
  while(!feof(fp_in)) {
    if(fgets(in_line, 80, fp_in)) {
      if(in_line[0] != '%') break;  
    }
  }

/* Read orbital elements from input file: */
/* Marco's format, i.e., orbit_format = -1 or orbit_format = 1 */
   if(orbit_format != -1 && orbit_format != 1) {
      fprintf(stderr, "read_orbital_elements_from_file/Fatal error: bad option: orbit_format = %d\n", 
              orbit_format);
      exit(-1);
      }
/* JLP2010: I add the possibility of errors */
   if(!strncmp(in_line, "#Errors: ",9)) {
   expected_nval = 7;
   nval = sscanf(in_line, "#Errors: %lf %lf %lf %lf %lf %lf %lf", 
                 &err_Omega_node[k], &err_omega_peri[k], &err_i_incl[k], 
                 &err_e_eccent[k], &err_T_periastron[k], &err_Period[k], 
                 &err_a_smaxis[k]);
   if(nval != expected_nval) {
     fprintf(stderr, "read_orbital_elements_from_file/Error reading errors from input file (nval=%d expected_nval=%d) \n %s \n", 
              nval, expected_nval, in_line);
     fclose(fp_in);
     return(1);
    }

#ifdef DEBUG
 printf("nval=%d err_Omega_node=%.3f err_omega_peri=%.3f err_incl=%.3f err_e=%.3f err_T=%.3f err_P=%.3f err_a=%.3f\n", 
        nval, err_Omega_node[k], err_omega_peri[k], err_i_incl[k], 
        err_e_eccent[k], err_T_periastron[k], err_Period[k], err_a_smaxis[k]);
#endif
/* Conversion to radians (necessary for ephemerids): */
  err_Omega_node[k] *= DEGTORAD;
  err_omega_peri[k] *= DEGTORAD;
  err_i_incl[k] *= DEGTORAD;

/***************************************************************************/
   } else {
   k++;
   expected_nval = 8;
   nval = sscanf(in_line, "%lf %lf %lf %lf %lf %lf %lf %lf", 
                 &Omega_node[k], &omega_peri[k], &i_incl[k], &e_eccent[k], 
                 &T_periastron[k], &Period[k], &a_smaxis[k], &orbit_equinox[k]);
   if(nval != expected_nval) {
     fprintf(stderr, "read_orbital_elements_from_file/Fatal error reading input file (nval=%d expected_nval=%d) \n %s \n", 
              nval, expected_nval, in_line);
     fclose(fp_in);
     return(-2);
    }

#ifdef DEBUG
 printf("nval=%d Omega_node=%.3f omega_peri=%.3f incl=%.3f e=%.3f T=%.3f P=%.3f a=%.3f Equinox=%.3f\n", 
        nval, Omega_node[k], omega_peri[k], i_incl[k], e_eccent[k], 
        T_periastron[k], Period[k], a_smaxis[k], orbit_equinox[k]);
#endif

/* Conversion to radians (necessary for ephemerids): */
  Omega_node[k] *= DEGTORAD;
  omega_peri[k] *= DEGTORAD;
  i_incl[k] *= DEGTORAD;
  mean_motion[k] = (360.0 / Period[k]) * DEGTORAD;
  }
} /* EOF loop on k */


/* Compute error on mean motion: */
for(k = 0; k < nber_of_orbits; k++) {
  err_mean_motion[k] = mean_motion[k] * err_Period[k] / Period[k];
  }

 fclose(fp_in);

if(k == nber_of_orbits) 
 status = 0;
else
 status = -1;

return(status);
}
/***************************************************************************
* Compute Thiele elements from orbital elements
*
* Omega_node (radians)
* omega_peri (radians)
* i_incl (radians)
* a_smaxis (arcsec): semi-major-axis
****************************************************************************/
int compute_Thiele_elements(double Omega_node, double omega_peri, 
                            double i_incl, double e_eccent, double T_periastron,
                            double Period, double a_smaxis, double mean_motion,
                            double orbit_equinox, double *AA, double *BB, 
                            double *FF, double *GG)
{
*AA = a_smaxis * (cos(omega_peri) * cos(Omega_node) 
                - sin(omega_peri) * sin(Omega_node) * cos(i_incl));
*BB = a_smaxis * (cos(omega_peri) * sin(Omega_node) 
                + sin(omega_peri) * cos(Omega_node) * cos(i_incl));
*FF = a_smaxis * (-sin(omega_peri) * cos(Omega_node) 
                - cos(omega_peri) * sin(Omega_node) * cos(i_incl));
*GG = a_smaxis * (-sin(omega_peri) * sin(Omega_node) 
                + cos(omega_peri) * cos(Omega_node) * cos(i_incl));
return(0);
}
