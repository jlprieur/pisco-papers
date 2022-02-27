/************************************************************************
* "PISCO_catalog_utils.h"
* Set of routines used to read PISCO catalog
*
* JLP 
* Version 21/06/2015
*************************************************************************/
#ifndef _PISCO_catalog_utils_h /* BOF sentry */
#define _PISCO_catalog_utils_h

#include <stdio.h>
#include <stdlib.h>                  /* exit() */
#include <string.h>
#include <ctype.h>                   /* isprint... */
#include <math.h>
#include <time.h>                    /* date */

#ifndef MAXI
#define MAXI(a,b) ((a) < (b)) ? (b) : (a)
#endif

#ifndef ABS
#define ABS(a) ((a) < 0.0  ? (-(a)) : (a))
#endif

#ifndef PI
#define PI 3.14159265
#endif

#define DEGTORAD   (PI/180.00)

/* Declaring linkage specification to have "correct names"
* that can be linked with C programs */

#ifdef __cplusplus
extern "C" {
#endif

/* PISCO catalog zeiss_doppie.cat with the list of objects used by TAV1.EXE */
int get_coordinates_from_PISCO_catalog(char *PISCO_catalog_name, 
                                       char *NameInPiscoCatalog,
                                       double *alpha, double *delta, 
                                       double *coord_equinox); 
int get_data_from_PISCO_catalog(char *PISCO_catalog_name, 
                                char *NameInPiscoCatalog,
                                double *magV, double *B_V, 
                                double *paral, double *err_paral, 
                                double *magV_A, double *magV_B,
                                char *spectral_type, char *discov_name,
                                char *comp_name, char *ads_name, 
                                char *WDS_name);
int read_coordinates_from_PISCO_catalog(char *line_buffer, double *alpha, 
                                        double *delta, double *coord_equinox);
int PISCO_catalog_read_line1(char *in_line, char *object_name, char *comp_name,
                             char *ADS_name, double *alpha_Pcat, 
                             double *delta_Pcat, double *equinox_Pcat);
int PISCO_catalog_read_line2(char *in_line, double *magV_A, double *magV_B, 
                             char *spectral_type);
int PISCO_catalog_read_line2_discov(char *in_line, char *discov_name,
                                    char *comp_name);
int PISCO_catalog_read_line2_WDS(char *in_line, char *WDS_name);
int PISCO_catalog_read_line2_Hip(char *in_line, double *magV, double *B_V,     
                                 double *paral, double *err_paral);

#ifdef __cplusplus
}
#endif


#endif
