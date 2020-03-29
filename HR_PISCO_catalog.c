/************************************************************************
* "HR_PISCO_catalog.c"
* To create a Hertzprung-Russel diagram of the objects contained 
* in the PISCO catalog with the list of targets " zeiss_doppie_new.cat"
* (that is used by TAV1.EXE)
*
* JLP 
* Version 25/05/2009
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>                  /* exit() */
#include <string.h>
#include <ctype.h>                   /* isprint... */
#include <math.h>
#include <time.h>                    /* date */
#include "jlp_catalog_utils.h"       /* Routines used to read catalogs */

#define ABS(a) ((a) < 0.0  ? (-(a)) : (a))
#ifndef PI
#define PI 3.14159265
#endif
#define DEGTORAD   (PI/180.00)

/*
#define DEBUG
#define DEBUG_1
*/

static int create_HR_curve(FILE *fp_PISCO_cat, FILE *fp_out,
                           double paral_rel_error0);

int main(int argc, char *argv[])
{
char PISCO_catalog[80], HR_filename[80]; 
FILE *fp_PISCO_cat, *fp_out; 
time_t ttime = time(NULL);
double paral_rel_error0;

if(argc != 3 && argc != 4) {
  printf("argc = %d\n", argc);
  printf("Syntax: HR_PISCO_catalog PISCO_catalog HR_or_photom_file [relative_error_of_parallax]\n");
  printf("Example: HR_PISCO_catalog zeiss_doppie_new.cat HR1.dat 0.5\n");
  printf("Example: HR_PISCO_catalog zeiss_doppie_new.cat Photom1.dat -1.\n");
  printf("Enter a negative relative error to output photometric data of all the objects\n");
  printf("Enter a positive relative error to output absolute magnitudes of objects measured by Hipparcos\n");
  return(-1);
}
strcpy(PISCO_catalog, argv[1]);
strcpy(HR_filename, argv[2]);
if(argc == 4) {
 sscanf(argv[3], "%lf", &paral_rel_error0);
 } else {
 paral_rel_error0 = 0.8;
 }

printf("OK: cat=%s HR_file=%s \n", PISCO_catalog, HR_filename);

/* Open input file containing the Latex PISCO catalog */
if((fp_PISCO_cat = fopen(PISCO_catalog, "r")) == NULL) {
   fprintf(stderr, "HR_PISCO_catalog/Fatal error opening PISCO catalog: %s\n",
           PISCO_catalog);
   return(-1);
  }

/* Open output file containing the HR diagram */
if((fp_out = fopen(HR_filename, "w")) == NULL) {
   fprintf(stderr, "HR_PISCO_catalog/Fatal error opening output HR file: %s\n",
           HR_filename);
   fclose(fp_PISCO_cat);
   return(-1);
  }

/* Header of the output photometric file: */
if(paral_rel_error0 < 0) {
  fprintf(fp_out, "%% Photometric data of %s \n%% Created on %s", 
        PISCO_catalog, ctime(&ttime));
  fprintf(fp_out, "%% magV magV_A B-V Delta_mag discov_name\n");
/* Header of the output HR file: */
  } else {
  fprintf(fp_out, "%% HR diagram from %s (paral_rel_error0=%.3f) \n%% Created on %s", 
        PISCO_catalog, paral_rel_error0, ctime(&ttime));
  fprintf(fp_out, "%% magV magV_A B-V Delta_mag MagV MagV_A discov_name\n");
  }

/* Scan the PISCO catalog, 
* compute M_V, B-V, and build a HR diagram if paral_rel_error0 > 0
*/
create_HR_curve(fp_PISCO_cat, fp_out, paral_rel_error0); 

/* Close opened files:
*/
fclose(fp_PISCO_cat);
fclose(fp_out);
return(0);
}
/***************************************************************************
* Scan the PISCO catalog
* Example:
* ADS 326AB (3.40, 340) &  0 23 59.8 &-03 28 31 & 2000.0
* & A: 7.4, B: 9.3 A2 & BU   488AB, & WDS0000+000 HIP234 PI_HIP=0.34(+/-0.02) V_HIP=7.6 B-V_HIP=0.34 \cr
*
* INPUT:
* fp_PISCO_cat: pointer to the Latex PISCO catalog
* fp_out: pointer to the output HR file 
* paral_rel_error0: relative error on the parallax
*                   (if negative: only photometry is required)
*
* OUTPUT:
* 
***************************************************************************/
static int create_HR_curve(FILE *fp_PISCO_cat, FILE *fp_out,
                           double paral_rel_error0)
{
char in_line[128], spectral_type[40];
char object_name[40], comp_name[20], ADS_name[20], discov_name[20];
double alpha_Pcat, delta_Pcat, equinox_Pcat;
double magV_A, magV_B, magV, B_V, paral, err_paral, magV_abs, magV_A_abs;
int status, first_line_is_read, iline, n_objects, n_focusing;
int has_an_orbit;

/* Scan the input Latex PISCO catalog: */
first_line_is_read = 0;
iline = 0;
n_objects = 0;
/* Number of single stars used for focusing the telescope: */
n_focusing = 0;

/* DEBUG: */
#ifdef DEBUG
printf("DEBUG/WARNING: limitation to the first 100 lines! \n");
while(!feof(fp_PISCO_cat) && iline < 100) {
printf("DDEBUG: limited to 100 lines! \n");
#else
while(!feof(fp_PISCO_cat)) {
#endif
  if(fgets(in_line, 120, fp_PISCO_cat)) {
    iline++;
    if(in_line[0] != '%') {
/* Check if it is a first line or a second line */
      if(in_line[0] != '&') {
/* Read the first line */
         status = PISCO_catalog_read_line1(in_line, object_name, comp_name, 
                                           ADS_name, &alpha_Pcat, &delta_Pcat, 
                                           &equinox_Pcat);
/* Case of single stars used during observations for focusing the telescope: */
         if(status == 10) {
           n_focusing++;
           break;
           }

         if(status) {
           fprintf(stderr, "Fatal error/Bad object syntax in line #%d\n", iline);
           exit(-1);
           }
#ifdef DEBUG_1
if(comp_name[0] != '\0') printf("iline=%d object name=>%s< comp_name=>%s<\n", 
                iline, object_name, comp_name);
#endif
         first_line_is_read = 1;
/* Read the second line */
      } else {
         if(!first_line_is_read) {
            fprintf(stderr, "Fatal error/Bad syntax in line #%d\n", iline);
            exit(-1);
            }
         first_line_is_read = 0;
         PISCO_catalog_read_line2_discov(in_line, discov_name);
         status = PISCO_catalog_read_line2(in_line, &magV_A, &magV_B,
                                           spectral_type, &has_an_orbit);
         if (status) {
            fprintf(stderr, "Fatal error/Bad syntax in line #%d\n", iline);
            exit(-1);
            } 
 
            status = PISCO_catalog_read_line2_Hip(in_line, &magV, &B_V,
                                                   &paral, &err_paral);
            if (status) {
              fprintf(stderr, "Fatal error/Bad syntax in line #%d\n", iline);
              exit(-1);
              } 
            if(paral > 0. && magV != 100.  && magV_A != 100.  
               && magV_B != 100. && B_V != 100.) {
/* Photometric data: */
               if(paral_rel_error0 < 0) {
                 compact_string(discov_name, 20);
/* % magV magV_A B-V Delta_mag discov_name spectral type*/
                 fprintf(fp_out, "%.3f %.3f %.3f %.3f %s %s\n", 
                         magV, magV_A, B_V, ABS(magV_A - magV_B), 
                         discov_name, spectral_type);
/* HR diagram: */
               } else if(err_paral/paral < paral_rel_error0) {
/* MV = mV - 2.5 log10 (d/d_10)^2
*  MV = mV - 5 log10 (d) + 5 log10(d_10)
*  MV = mV - 5 log10(d) + 5
* paral = 1/d
* MV = mV + 5 log10(paral) + 5
*/
/* Parallax in mas: */
                 magV_abs = magV + 5. * log10(paral * 0.001) + 5.;
                 magV_A_abs = magV_A + 5. * log10(paral * 0.001) + 5.;
                 compact_string(discov_name, 20);
/* % magV magV_A B-V Delta_mag MagV MagV_A discov_name spectral type*/
                 fprintf(fp_out, "%.3f %.3f %.3f %.3f %.3f %.3f %s %s\n", 
                         magV, magV_A, B_V, ABS(magV_A - magV_B), magV_abs, 
                         magV_A_abs, discov_name, spectral_type);
                 } /* else if(err_paral/paral < paral_rel_error0) */
              } /* EOF paral > 0 */
            n_objects++;
#ifdef DEBUG_1
if(comp_name[0] != '\0') printf("iline=%d discov name=>%s<\n", iline, 
                                 discov_name);
#endif
/* & A: 7.4, B: 9.3 A2 & BU   488AB, & WDS0000+000 HIP234 PI_HIP=0.34(+/-0.02) V_HIP=7.6 B-V_HIP=0.34 \cr
*/
      } /* EOF second line, starting with "&" */
    } /* EOF inline != %*/
  } /* EOF fgets */
} /* EOF while */

printf("End of process: %d lines read (%d confirmed binaries)\n", 
        iline, n_objects);
if(n_focusing) printf("... Exit at line #%d when focusing star was found\n",
                      iline);

return(0);
}
