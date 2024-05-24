/************************************************************************
* "residuals_2.c"
*
* To compute the residuals of all the measurements of a given binary star
* in order to evaluate an orbit.
*
* From omenc.for (version of 2008)
*
* JLP 
* Version 07/08/2009
*************************************************************************/
#include <math.h>
#include "jlp_catalog_utils.h"
#include "jlp_string.h"  // jlp_trim_string
#include "residuals_utils.h"

#define SQUARE(a) ((a)*(a))
#define DEBUG0
/*
#define DEBUG
*/
#define DEGTORAD   (PI/180.00)

static int scan_input_file(FILE *fp_in, FILE *fp_Drho_Dtheta, 
              FILE *fp_latex, FILE *fp_Dx_Dy, FILE *fp_orbit_data, 
              int extended_input_format, double *Omega_node, double *omega_peri,
              double *i_incl, double *e_eccent, double *T_periastron, 
              double *Period, double *a_smaxis, double *mean_motion, 
              double *orbit_equinox, int nber_of_orbits, double *sigma_rho,
              double *sigma_theta);
static int remove_bad_values(FILE *fp_in, FILE *fp_selection,
              int extended_input_format, double *Omega_node, double *omega_peri,
              double *i_incl, double *e_eccent, double *T_periastron, 
              double *Period, double *a_smaxis, double *mean_motion, 
              double *orbit_equinox, int nber_of_orbits, double sigma_rho_max,
              double sigma_theta_max);
static int residuals2_main(char* measures_infile, char *output_ext, 
              int extended_input_format, double *Omega_node, double *omega_peri,
              double *i_incl, double *e_eccent, double *T_periastron, 
              double *Period, double *a_smaxis, double *mean_motion, 
              double *orbit_equinox, int nber_of_orbits, double sigma_rho_max,
              double sigma_theta_max);

int main(int argc, char *argv[])
{
double Omega_node[3], omega_peri[3], i_incl[3], e_eccent[3]; 
double T_periastron[3], Period[3], a_smaxis[3], mean_motion[3], orbit_equinox[3];
double sigma_rho_max, sigma_theta_max;
char measures_infile[124], orbit_infile[124],output_ext[40];
int orbit_format, extended_input_format, nber_of_orbits;
int nval, status;

if(argc == 7) {
  if(*argv[6]) argc = 7;
  else if(*argv[5]) argc = 6;
  else if(*argv[4]) argc = 5;
  else if(*argv[3]) argc = 4;
  else if(*argv[2]) argc = 3;
  else if(*argv[1]) argc = 2;
  else argc = 1;
}
if(argc != 5 && argc != 6) {
  printf("Syntax: residuals_2 list_of_measures file_with_orbit orbit_format output_extension [nber of orbits,sigma_rho_max,sigma_theta_max]\n");
  printf("Format: 1 if Marco's format (Omega=node, omep=longitude of periastron, i, e, T, P, a, [equinox])\n");
  printf("        2 if OC6 format (P, a, i, Omega, T, e, omep, [equinox] )\n");
  printf("\n NB: sigma_rho_max and sigma_theta_max are only used for building the *_select file\n");
  printf("\n Enter sigma_rho_max=sigma_theta_max=-1 for selecting F,G,H,S codes (photo. and speckle)\n");
  return(-1);
}
strcpy(measures_infile, argv[1]);
strcpy(orbit_infile, argv[2]);
sscanf(argv[3], "%d", &orbit_format);
strcpy(output_ext, argv[4]);
if(argc == 6) {
 nval = sscanf(argv[5], "%d,%lf,%lf", &nber_of_orbits, &sigma_rho_max, 
               &sigma_theta_max);
 if(nval != 3 ) {
   fprintf(stderr, "Fatal error: bad syntax!\n");
   return(-1);
   }
 if(nber_of_orbits < 1 || nber_of_orbits > 3) {
   fprintf(stderr, "Fatal error: nber_of_orbits = %d not allowed here!\n",
           nber_of_orbits);
   return(-1);
   }
 } else {
 nber_of_orbits = 1;
 sigma_rho_max = 0.;
 sigma_theta_max = 0.;
 }

#ifdef DEBUG
printf("OK: measures=%s orbit=%s orbit_format=%d\n",
       measures_infile, orbit_infile, orbit_format);
printf("    output_ext=%s nber_of_orbits=%d\n", 
       output_ext, nber_of_orbits);
#endif

/* Read orbital parameters: 
*/
status = read_orbital_elements_from_file(orbit_infile, orbit_format, 
                                 Omega_node, omega_peri, i_incl, e_eccent, 
                                 T_periastron, Period, a_smaxis, mean_motion, 
                                 orbit_equinox, nber_of_orbits); 
if(status) {
 fprintf(stderr, "Fatal error reading orbital parameters!\n");
 return(-1);
 }
printf("JLPPPP: Omega=%f\n", Omega_node[0]);

/* Call residuals2_main that does the main job: 
*/
extended_input_format = 1;
residuals2_main(measures_infile, output_ext, extended_input_format,
                Omega_node, omega_peri, i_incl, e_eccent, T_periastron, Period,
                a_smaxis, mean_motion, orbit_equinox, nber_of_orbits, 
                sigma_rho_max, sigma_theta_max);

return(0);
}
/************************************************************************
* residuals2_main
* main routine of "residuals_2.c"
*
* INPUT:
* measures_infile: name of the file containing the (corrected) measurements
*                  as created by 1bin.for or orbit_weight.c
* output_ext: extension of the output files
* orbital elements of the orbit to be evaluated 
*
*************************************************************************/
static int residuals2_main(char* measures_infile, char *output_ext, 
              int extended_input_format, double *Omega_node, double *omega_peri,
              double *i_incl, double *e_eccent, double *T_periastron, 
              double *Period, double *a_smaxis, double *mean_motion, 
              double *orbit_equinox, int nber_of_orbits, double sigma_rho_max,
              double sigma_theta_max) 
{
double AA, BB, FF, GG, nn[3], sigma_rho, sigma_theta;
int k;
char filename[80];
FILE *fp_in, *fp_Drho_Dtheta, *fp_latex, *fp_Dx_Dy, *fp_selection;
FILE *fp_orbit_data;

fp_in = NULL;
fp_Drho_Dtheta = NULL;
fp_latex = NULL;
fp_Dx_Dy = NULL;
fp_selection = NULL;
fp_orbit_data = NULL;

// Open output data file containing the measurements, the residual vectors 
sprintf(filename, "%s_orb_resi.dat", output_ext);
jlp_trim_string(filename, 80);

// and the orbit (in XY plane) 
if((fp_orbit_data = fopen(filename, "w")) == NULL) {
   fprintf(stderr, "residuals2_main/Fatal error opening output file: %s\n",
           filename);
    return(-1);
  }
fprintf(fp_orbit_data, "%% residuals_2 (orbit_data) -- version 21/07/2018\n");
fprintf(fp_orbit_data, "%% x_O, y_O, x_C, y_C\n");


/* Open input file containing the measurements and the orbital parameters: */
if((fp_in = fopen(measures_infile, "r")) == NULL) {
   fprintf(stderr, "residuals2_main/Fatal error opening input file: %s\n",
           measures_infile);
    return(-1);
  }

/* Open output file with the measurements and the O-C residuals: */
sprintf(filename, "%s.OM", output_ext);
jlp_trim_string(filename, 80);

if((fp_Drho_Dtheta = fopen(filename, "w")) == NULL) {
   fprintf(stderr, "residuals2_main/Fatal error opening output file: %s\n",
           filename);
    fclose(fp_in);
    return(-1);
  }

/* Open output file with the measurements and the O-C residuals: */
sprintf(filename, "%s.dxy", output_ext);
jlp_trim_string(filename, 80);

/* Open file with residuals (Dx, Dy) to test whether there is a third body 
*/
 if((fp_Dx_Dy = fopen(filename, "w")) == NULL) {
   fprintf(stderr, "residuals2_main/Fatal error opening output file: %s\n",
           filename);
    if(fp_Drho_Dtheta != NULL) fclose(fp_Drho_Dtheta);
    fclose(fp_in);
    return(-1);
  }

/* Open latex output table with the residuals (and authors): */
sprintf(filename, "%s_resi.tex", output_ext);
jlp_trim_string(filename, 80);

if((fp_latex = fopen(filename, "w")) == NULL) {
   fprintf(stderr, "residuals2_main/Fatal error opening output file: %s\n",
           filename);
    if(fp_Dx_Dy != NULL) fclose(fp_Dx_Dy);
    if(fp_Drho_Dtheta != NULL) fclose(fp_Drho_Dtheta);
    fclose(fp_in);
    return(-1);
  }

for(k = 0; k < nber_of_orbits; k++) nn[k] = 360.0 / Period[k];

fprintf(fp_Drho_Dtheta, "%% residuals_2 (rho,theta) -- version 07/08/2009\n");
fprintf(fp_Dx_Dy, "%% residuals_2 (Dx,Dy) -- version 07/08/2009\n");
fprintf(fp_Drho_Dtheta, "%% Input file: %s\n", measures_infile);
fprintf(fp_Dx_Dy, "%% Input file: %s\n", measures_infile);

for(k = 0; k < nber_of_orbits; k++) {
fprintf(fp_Drho_Dtheta, "%% Orbital elements:\n");
fprintf(fp_Dx_Dy, "%% Orbital elements:\n");
fprintf(fp_Drho_Dtheta, "%% Omega_node=%.3f omega_peri=%.3f incl=%.3f e=%.4f \n", 
        Omega_node[k]/DEGTORAD, omega_peri[k]/DEGTORAD, i_incl[k]/DEGTORAD, 
        e_eccent[k]);
fprintf(fp_Dx_Dy, "%% Omega_node=%.3f omega_peri=%.3f incl=%.3f e=%.4f \n", 
        Omega_node[k]/DEGTORAD, omega_peri[k]/DEGTORAD, i_incl[k]/DEGTORAD, 
        e_eccent[k]);
fprintf(fp_Drho_Dtheta, "%% T=%.3f P=%.3f n=%.5f a=%.3f Equinox=%.3f\n", 
        T_periastron[k], Period[k], nn[k], a_smaxis[k], orbit_equinox[k]);
fprintf(fp_Dx_Dy, "%% T=%.3f P=%.3f n=%.5f a=%.3f Equinox=%.3f\n", 
        T_periastron[k], Period[k], nn[k], a_smaxis[k], orbit_equinox[k]);
compute_Thiele_elements(Omega_node[k], omega_peri[k], i_incl[k], e_eccent[k], 
                        T_periastron[k], Period[k], a_smaxis[k], mean_motion[k],
                        orbit_equinox[k], &AA, &BB, &FF, &GG);
fprintf(fp_Drho_Dtheta, "%% Thiele elements: A=%12.5f B=%12.5f F=%12.5f G=%12.5f\n",        AA, BB, FF, GG); 
} /* EOF loop on k */

if(extended_input_format)
  fprintf(fp_Dx_Dy, "%% epoch Dx Dy nights author aperture weight\n");
else
  fprintf(fp_Dx_Dy, "%% epoch Dx Dy weight\n");

fprintf(fp_Drho_Dtheta, "%% epoch rho_O rho_C Drho_O-C theta_O theta_C Dtheta_O-C author\n");

fprintf(fp_latex, "\\begin{tabular}{crrl}\n");
fprintf(fp_latex, "\\hline\n");
fprintf(fp_latex, "Epoch & $\\Delta \\rho$ (O-C) & $\\Delta \\theta$ (O-C) & Observer \\\\\n");
fprintf(fp_latex, " &  (\\arcsec) & (\\degr) & \\\\\n");
fprintf(fp_latex, "\\hline\n");

/* Scan the input file, compute the residuals and write results to output files
*/
scan_input_file(fp_in, fp_Drho_Dtheta, fp_latex, fp_Dx_Dy, fp_orbit_data, 
                extended_input_format, Omega_node, omega_peri, i_incl, 
                e_eccent, T_periastron, Period, a_smaxis, mean_motion, 
                orbit_equinox, nber_of_orbits, &sigma_rho, &sigma_theta);
fclose(fp_in);
fclose(fp_Drho_Dtheta);
fclose(fp_latex);
fclose(fp_Dx_Dy);
fclose(fp_orbit_data);

/* Neutralize bad values: */
/* Re-open input file containing the measurements and the orbital parameters: */
if(sigma_rho_max != 0. || sigma_theta_max != 0.) {
if((fp_in = fopen(measures_infile, "r")) == NULL) {
   fprintf(stderr, "residuals2_main/Fatal error re-opening input file: %s\n",
           measures_infile);
    return(-1);
  }

/* Open output file which is a copy of the input file
* with bad measurements removed : */
sprintf(filename, "%s_select.dat", output_ext);
jlp_trim_string(filename, 80);

 if((fp_selection = fopen(filename, "w")) == NULL) {
   fprintf(stderr, "residuals2_main/Fatal error opening output file: %s\n",
           filename);
    fclose(fp_in);
    return(-1);
  }

remove_bad_values(fp_in, fp_selection, 
                extended_input_format, Omega_node, omega_peri, i_incl, 
                e_eccent, T_periastron, Period, a_smaxis, mean_motion, 
                orbit_equinox, nber_of_orbits, sigma_rho_max, 
                sigma_theta_max);

fclose(fp_in);
fclose(fp_selection);
} else {
 printf("sigma_rho_max = %f sigma_theta_max = %f, hence no selection!\n",
         sigma_rho_max, sigma_theta_max);
}


return(0);
}
/********************************************************************
* Scan the input file, compute the residuals and write results to output files
*
* OUTPUT:
*  sigma_rho, sigma_theta: standard deviation in rho and theta of the residuals
*********************************************************************/
static int scan_input_file(FILE *fp_in, FILE *fp_Drho_Dtheta,  
              FILE *fp_latex, FILE *fp_Dx_Dy, FILE *fp_orbit_data, 
              int extended_input_format, double *Omega_node, double *omega_peri,
              double *i_incl, double *e_eccent, double *T_periastron, 
              double *Period, double *a_smaxis, double *mean_motion, 
              double *orbit_equinox, int nber_of_orbits, double *sigma_rho,
              double *sigma_theta) 
{
double epoch, rho_o, theta_o, weight, rho_c, theta_c, Drho, Dtheta; 
double c_tolerance[3], Drho_err, Dtheta_err, x_O, y_O, x_C, y_C;
double Drho_sum, Drho_sumsq, Dtheta_sum, Dtheta_sumsq, Dx, Dy;
int n_nights, iaperture, n_observations, iline, nval, n_Drho, n_Dtheta, k;
char buffer[80], author[10];

/*  c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*/
for(k = 0; k < nber_of_orbits; k++) 
c_tolerance[k] = ABS(1.5E-5 * cos(i_incl[k]) 
                   / sqrt((1.0 + e_eccent[k])/(1.0 - e_eccent[k])));

n_observations = 0;
iline = 0;
Drho_sum = 0.;
Drho_sumsq = 0.;
Dtheta_sum = 0.;
Dtheta_sumsq = 0.;
n_Drho = 0;
n_Dtheta = 0;

/* Main loop: */
while(!feof(fp_in)) {
  if(!fgets(buffer, 80, fp_in)) break;
  iline++;
/* Possibility of commented lines, starting with % or # : */
  if(buffer[0] != '%' && buffer[0] != '#') {
    if(extended_input_format) {
      nval = sscanf(buffer, " %8lf %8lf %8lf %2d %3s %d %lf\n",
             &epoch, &rho_o, &theta_o, &n_nights, author, &iaperture, &weight);
      if(nval < 7) {
       fprintf(stderr, "Error reading line %d\n", iline);
       break;
       }
/* Weight is already present in non-extended format input file: */
     } else {
      nval = sscanf(buffer, " %8lf %8lf %8lf %lf\n", &epoch, &rho_o, &theta_o, 
                     &weight);
      n_nights = 0; author[0] = '\0'; iaperture = 0;
      if(nval != 4) {
        fprintf(stderr, "Error reading line %d\n", iline);
        break;
        }
     }
/* Exit from loop (and exit from program) when Epoch = 0.0 is found:
*/
      if(epoch == 0.0) break;

#ifdef DEBUG
 if(n_observations < 10)
    printf("epoch=%.2f rho_o=%.3f theta_o=%.1f author=%3s weight=%.2f\n",
            epoch, rho_o, theta_o, author, weight);
#endif

/* Compute the ephemerids corresponding to the observation epoch: */
compute_ephemerid_of_multiple_system(nber_of_orbits,
                      Omega_node, omega_peri, i_incl, e_eccent, T_periastron, 
                      Period, a_smaxis, mean_motion, epoch, c_tolerance, 
                      &theta_c, &rho_c);

 Dtheta = theta_o - theta_c;
 if(Dtheta < -300.0) Dtheta += 360.0;
 if(Dtheta > 300.0) Dtheta -= 360.0;
 n_Dtheta++;
 Dtheta_sum += Dtheta;
 Dtheta_sumsq += Dtheta*Dtheta;

 if(rho_o > 0) {
   Drho = rho_o - rho_c;
   Drho_sum += Drho;
   Drho_sumsq += Drho*Drho;
   n_Drho++;
   } else {
   rho_o = -100.;
   Drho = -100.;
   }
/* epoch rho_O rho_C Drho_O-C theta_O theta_C Dtheta_O-C author
*/
 fprintf(fp_Drho_Dtheta, "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %-3s\n",
         epoch, rho_o, rho_c, Drho, theta_o, theta_c, Dtheta, author);
/* Use $ $ in order to have long minus signs with LaTeX... */
 if(Drho == -100) {
  fprintf(fp_latex, "%8.3f & \\nodata & $%.3f$ & %s \\\\\n",
          epoch, Dtheta, author);
  } else {
  fprintf(fp_latex, "%8.3f & $%.3f$ & $%.3f$ & %s \\\\\n",
          epoch, Drho, Dtheta, author);
  }

/* Output the residuals as (Dx,Dy) */
 if(Drho != -100) {
/* DEGTORAD = (PI/180.00) */
    x_O = rho_o * cos(theta_o * DEGTORAD);
    x_C = rho_c * cos(theta_c * DEGTORAD);
    Dx = x_O - x_C;
    y_O = rho_o * sin(theta_o * DEGTORAD);
    y_C = rho_c * sin(theta_c * DEGTORAD);
    Dy = y_O - y_C;
// orbit_data : x_O, y_O, x_C, y_C
    fprintf(fp_orbit_data, " %8.3f %8.3f %8.3f %8.3f\n",
            x_O, y_O, x_C, y_C);
    if(extended_input_format) 
      fprintf(fp_Dx_Dy, " %8.3f %8.3f %8.3f %2d %-3s %3d %4.1f\n",
              epoch, Dx, Dy, n_nights, author, iaperture, weight);
    else 
      fprintf(fp_Dx_Dy, " %8.3f %8.3f %8.3f %4.1f\n", epoch, Dx, Dy, 
              weight);
  }

 n_observations++;
 } /* EOF not commented line */

} /* EOF while feof */

Drho_err = 0.;
Dtheta_err = 0.;
if(n_Drho > 3) {
  Drho_sum /= (double)n_Drho;
  Drho_sumsq /= (double)n_Drho;
  Drho_err = sqrt(Drho_sumsq - Drho_sum * Drho_sum);
  }
if(n_Dtheta > 3) {
  Dtheta_sum /= (double)n_Dtheta;
  Dtheta_sumsq /= (double)n_Dtheta;
  Dtheta_err = sqrt(Dtheta_sumsq - Dtheta_sum * Dtheta_sum);
  }
printf("n_observations =%d (n_Dtheta=%d n_Drho=%d)\n", 
        n_observations, n_Dtheta, n_Drho);
printf("mean error rms: Drho_O-C=%.3f Dtheta_O-C=%.3f\n", 
        Drho_err, Dtheta_err);

/* Store results to output variables: */
*sigma_rho = Drho_err;
*sigma_theta = Dtheta_err;

fprintf(fp_Drho_Dtheta, "%% n_observations =%d (n_Dtheta=%d n_Drho=%d)\n", 
        n_observations, n_Dtheta, n_Drho);
fprintf(fp_Drho_Dtheta, "%% mean error rms: Drho_O-C=%.3f Dtheta_O-C=%.3f\n", 
        Drho_err, Dtheta_err);

fprintf(fp_latex, "\\hline\n");
fprintf(fp_latex, "\\end{tabular}\n");

return(0);
}
/********************************************************************
* Remove all aberrant values if residuals are above a given threshold
* sigma_rho_max, sigma_theta_max
*        if -1,-1 selection of codes F,G,H,S (photographic and speckle)
*********************************************************************/
static int remove_bad_values(FILE *fp_in, FILE *fp_selection,
              int extended_input_format, double *Omega_node, double *omega_peri,
              double *i_incl, double *e_eccent, double *T_periastron, 
              double *Period, double *a_smaxis, double *mean_motion, 
              double *orbit_equinox, int nber_of_orbits, double sigma_rho_max,
              double sigma_theta_max) 
{
double epoch, rho_o, theta_o, weight, rho_c, theta_c, Drho, Dtheta; 
double c_tolerance[3];
int n_nights, iaperture, n_full_observations, n_observations; 
int iline, nval, k, code_selection;
char buffer[80], author[10], code[2];
double sumsq_dtheta, sumsq_drho;

/*  if -1,-1 selection of codes F,G,H,S (photographic and speckle)
*/
code_selection = (sigma_rho_max == -1 || sigma_theta_max == -1) ? 1 : 0;

/*  c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*/
for(k = 0; k < nber_of_orbits; k++) 
c_tolerance[k] = ABS(1.5E-5 * cos(i_incl[k]) 
                   / sqrt((1.0 + e_eccent[k])/(1.0 - e_eccent[k])));

n_observations = 0;
sumsq_drho = 0.;
sumsq_dtheta = 0.;
n_full_observations = 0;
iline = 0;

/* Main loop: */
while(!feof(fp_in)) {
  if(!fgets(buffer, 80, fp_in)) break;
  iline++;
/* Simply copy input to output if commented line: */
  if(buffer[0] == '%' || buffer[0] == '#') {
  fprintf(fp_selection, buffer);
  } else {
    if(extended_input_format) {
      nval = sscanf(buffer, " %8f %8f %8f %2d %3s %d %f %c\n",
             &epoch, &rho_o, &theta_o, &n_nights, author, &iaperture, &weight,
             code);
      if(nval == 7) code[0] = ' ';
      if(nval < 7) {
       fprintf(stderr, "Error reading line %d\n", iline);
       break;
       }
/* Weight is already present in non-extended format input file: */
     } else {
      nval = sscanf(buffer, " %8f %8f %8f %f\n", &epoch, &rho_o, &theta_o, 
                     &weight);
      n_nights = 0; author[0] = '\0'; iaperture = 0; code[0] = ' ';
      if(nval != 4) {
        fprintf(stderr, "Error reading line %d\n", iline);
        break;
        }
     }
/* Exit from loop (and exit from program) when Epoch = 0.0 is found:
*/
   if(epoch == 0.0) break;

   code[1] = '\0';
#ifdef DEBUG0
 if(n_observations < 10) 
    printf("epoch=%.2f rho_o=%.3f theta_o=%.1f author=%3s weight=%.2f code=%s\n",
            epoch, rho_o, theta_o, author, weight, code);
#endif

if(code_selection) {
/* Copy input to output if code is OK: */
   switch (code[0]) { 
     case 'F':
     case 'G':
     case 'H':
     case 'S':
       fprintf(fp_selection, buffer);
       n_observations++;
       if(rho_o > 0 ) n_full_observations++;
       break;
     default:
       break;
     }
} else {
/* Compute the ephemerids corresponding to the observation epoch: */
compute_ephemerid_of_multiple_system(nber_of_orbits,
                      Omega_node, omega_peri, i_incl, e_eccent, T_periastron, 
                      Period, a_smaxis, mean_motion, epoch, c_tolerance, 
                      &theta_c, &rho_c);

 Drho = rho_o - rho_c;
/* When rho is negative, rho values are not taken into account 
* for orbit computation BUT the associated angular measurements are!)
*/
 if(rho_o < 0) Drho = 0.;
 Dtheta = theta_o - theta_c;
 if(Dtheta < -300.0) Dtheta += 360.0;
 if(Dtheta > 300.0) Dtheta -= 360.0;

#ifdef DEBUG0
 if(n_observations < 0)
    printf("Drho=%.3f Dtheta=%.1f (sigma_rho_max=%f theta_max=%f)\n", Drho, Dtheta,
           sigma_rho_max, sigma_theta_max);
#endif

/* Copy input to output if small residual: */
   if(SQUARE(Dtheta) <= SQUARE(sigma_theta_max) 
     && SQUARE(Drho) <= SQUARE(sigma_rho_max)){
      fprintf(fp_selection, buffer);
      n_observations++;
      sumsq_drho += SQUARE(Drho);
      sumsq_dtheta += SQUARE(Dtheta);
      if(rho_o > 0 ) n_full_observations++;
/* Neutralize rho by adding a minus sign 
* if small residual on theta and large residual in rho
*/
   } else if(SQUARE(Dtheta) <= SQUARE(sigma_theta_max) 
     && SQUARE(Drho) > SQUARE(sigma_rho_max)){

/* Extended format:
* author is limited to the first 3 characters 
* (and is left justified with %-3s) 
*/
     n_observations++;
     sumsq_drho += SQUARE(Drho);
     sumsq_dtheta += SQUARE(Dtheta);
     fprintf(fp_selection," %8.3f %8.3f %8.3f %2d %-3s %3d %4.1f %1s\n",
             epoch, -rho_o, theta_o, n_nights, author, iaperture, 
             weight, code);

/* DEBUG: */
    printf("GOOD THETA and BAD VALUE FOR RHO: epoch=%.2f rho_o=%.3f theta_o=%.1f author=%3s weight=%.2f code=%s\n",
            epoch, rho_o, theta_o, author, weight, code);
    printf("rho_c=%.3f theta_c=%.1f Drho=%.3f Dtheta=%.3f\n",
            rho_c, theta_c, Drho, Dtheta);
   } /* EOF case of bad value */
  } /* EOF case of sigma selection */
 } /* EOF not commented line */

} /* EOF while feof */

if(n_observations > 0) {
 printf("After selection: n_observations=%d mean Drho=%f Dtheta=%f\n",
      n_observations, sqrt(sumsq_drho / (double)n_observations),
      sqrt(sumsq_dtheta / (double)n_observations));
}

if(code_selection) {
printf("%% %d incomplete and %d complete observations selected with F,G,H,S code\n", 
       n_observations, n_full_observations);
fprintf(fp_selection,"%% %d incomplete and %d complete observations selected with F,G,H,S code\n", 
       n_observations, n_full_observations);
} else {
printf("%% %d incomplete and %d complete observations selected with residuals smaller than %.3f arcseconds and %.3f degrees\n", 
       n_observations, n_full_observations, sigma_rho_max, sigma_theta_max);
fprintf(fp_selection,"%% %d incomplete and %d complete observations selected with residuals smaller than %.3f degrees and %.3f arcseconds\n", 
       n_observations, n_full_observations, sigma_rho_max, sigma_theta_max);
}
return(0);
}
