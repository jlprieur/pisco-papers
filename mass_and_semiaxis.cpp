/************************************************************************
* "mass_and_semiaxis.c"
*
* To compute the mass and linear semi-major axis of a binary star
* whose orbit is known.
*
* JLP 
* Version 12/10/2010
*************************************************************************/
#include <stdio.h>
#include <math.h>
#include "jlp_catalog_utils.h"
#include "residuals_utils.h"

#define SQUARE(a) ((a)*(a))
#define DEBUG
#define DEGTORAD   (PI/180.00)

static int compute_mass_and_semiaxis(char *out_filename, double *Omega_node, 
              double *omega_peri, double *i_incl, double *e_eccent, 
              double *T_periastron, double *Period, double *a_smaxis, 
              double *mean_motion, double *orbit_equinox, double parallax, 
              double *err_Omega_node, double *err_omega_peri, double *err_i_incl,
              double *err_e_eccent, double *err_T_periastron, double *err_Period,
              double *err_a_smaxis, double *err_mean_motion, double err_parallax, 
              int nber_of_orbits);

int main(int argc, char *argv[])
{
double Omega_node[3], omega_peri[3], i_incl[3], e_eccent[3]; 
double T_periastron[3], Period[3], a_smaxis[3], mean_motion[3], orbit_equinox[3];
double err_Omega_node[3], err_omega_peri[3], err_i_incl[3], err_e_eccent[3]; 
double err_T_periastron[3], err_Period[3], err_a_smaxis[3], err_mean_motion[3]; 
double parallax, err_parallax;
int orbit_format, nber_of_orbits, nval, status;
char out_filename[80], orbit_infile[80];

if(argc == 7) {
  if(*argv[6]) argc = 7;
  else if(*argv[5]) argc = 6;
  else if(*argv[4]) argc = 5;
  else if(*argv[3]) argc = 4;
  else if(*argv[2]) argc = 3;
  else if(*argv[1]) argc = 2;
  else argc = 1;
}
if(argc != 5) {
  printf("Syntax: mass_and_semiaxis outfile infile_with_orbit orbit_format,number_of_orbits parallax_msec,err_parallax_msec\n");
  printf("Format: 1 if Marco's format (Omega=node, omep=longitude of periastron, i, e, T, P, a, [equinox])\n");
  printf("        2 if OC6 format (P, a, i, Omega, T, e, omep, [equinox] )\n");
  return(-1);
}
strcpy(out_filename, argv[1]);
strcpy(orbit_infile, argv[2]);
nval = sscanf(argv[3], "%d,%d", &orbit_format, &nber_of_orbits);
if(nval == 1) nber_of_orbits = 1;
sscanf(argv[4], "%lf,%lf", &parallax, &err_parallax);
if(nval == 1) err_parallax = 0.;

#ifdef DEBUG
printf("OK: orbit in %s orbit_format=%d\n",
       orbit_infile, orbit_format);
printf("    output to %s nber_of_orbits=%d\n", 
       out_filename, nber_of_orbits);
#endif

/* Read orbital parameters: 
*/
status = read_orbital_elements_and_errors_from_file(orbit_infile, orbit_format, 
                                 Omega_node, omega_peri, i_incl, e_eccent, 
                                 T_periastron, Period, a_smaxis, mean_motion, 
                                 orbit_equinox, err_Omega_node, err_omega_peri,
                                 err_i_incl, err_e_eccent, err_T_periastron, 
                                 err_Period, err_a_smaxis, err_mean_motion, 
                                 nber_of_orbits); 
if(status == 1) {
 fprintf(stderr, "Error reading errors of orbital parameters!\n");
 fprintf(stderr, "Will assume null errors.\n");
 } else if(status == -1) {
 fprintf(stderr, "Fatal error reading orbital parameters and their errors!\n");
 exit(-1);
 }

/* Call compute_mass_and_semiaxis that does the main job: 
*/
compute_mass_and_semiaxis(out_filename, Omega_node, omega_peri, i_incl, 
                          e_eccent, T_periastron, Period, a_smaxis, 
                          mean_motion, orbit_equinox, parallax,
                          err_Omega_node, err_omega_peri, err_i_incl, 
                          err_e_eccent, err_T_periastron, err_Period, 
                          err_a_smaxis, err_mean_motion, err_parallax, 
                          nber_of_orbits);

return(0);
}
/************************************************************************
* compute_mass_and_semiaxis 
* main routine of "residuals_2.c"
*
* INPUT:
* measures_infile: name of the file containing the (corrected) measurements
*                  as created by 1bin.for or orbit_weight.c
* output_ext: extension of the output files
* orbital elements of the orbit to be evaluated 
*
*************************************************************************/
static int compute_mass_and_semiaxis(char *out_filename, double *Omega_node, 
              double *omega_peri, double *i_incl, double *e_eccent, 
              double *T_periastron, double *Period, double *a_smaxis, 
              double *mean_motion, double *orbit_equinox, double parallax, 
              double *err_Omega_node, double *err_omega_peri, 
              double *err_i_incl, double *err_e_eccent, 
              double *err_T_periastron, double *err_Period,
              double *err_a_smaxis, double *err_mean_motion, 
              double err_parallax, int nber_of_orbits) 
{
double a_AU, err_a_AU, sum_masses, err_sum_masses;
FILE *fp_out;

/* Open output file containing the measurements and the orbital parameters: */
if((fp_out = fopen(out_filename, "w")) == NULL) {
   fprintf(stderr, "compute_mass_and_semiaxis/Fatal error opening output file: %s\n",
           out_filename);
    return(-1);
  }

printf("OK: Omega=%f err_Omega=%f\n", Omega_node[0], err_Omega_node[0]);

fprintf(fp_out, "Parallax = %.2f +/- %.2f (mas)\n", parallax, err_parallax); 

/* Semi major axis: conversion to astronomical units (AU) */
a_AU = a_smaxis[0] / (parallax / 1000.);
err_a_AU = a_AU * sqrt(SQUARE(err_a_smaxis[0] / a_smaxis[0]) 
                       + SQUARE(err_parallax / parallax));

printf("Semi-axis: a = %.2f +/- %.2f (Astron. Units)\n", 
        a_AU, err_a_AU);
fprintf(fp_out, "Semi-axis: a = %.2f +/- %.2f (Astron. Units)\n", 
         a_AU, err_a_AU);

/* Sum of the masses: computation of value and error */
sum_masses = pow(a_AU, 3) / (Period[0] * Period[0]);
err_sum_masses = sum_masses * sqrt( 9.* SQUARE(err_a_smaxis[0] / a_smaxis[0])
                               + 9.* SQUARE(err_parallax / parallax) 
                                + 4. * SQUARE(err_Period[0] / Period[0]));
printf("Sum of masses: M_1 + M_2 = %.2f +/- %.2f (Solar Masses)\n", 
        sum_masses, err_sum_masses);
fprintf(fp_out, "Sum of masses: M_1 + M_2 = %.2f +/- %.2f (Solar Masses)\n", 
        sum_masses, err_sum_masses);

fclose(fp_out);

return(0);
}
